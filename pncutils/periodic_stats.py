#!/usr/bin/env python3

import PseudoNetCDF as pnc
import numpy as np
import datetime
import warnings

known_spc = ['A24PM25', 'MDA8O3']

# fudge to work around bug for verdi...
#waste_first_frame = True
waste_first_frame = False

oned = datetime.timedelta(days=1)

def iodate_to_pydate(iodate):
    jd = iodate % 1000
    yr = iodate // 1000
    return datetime.date(yr, 1, 1) + oned * (jd-1)
def pydate_to_iodate(pydate):
    return pydate.year * 1000 + pydate.timetuple().tm_yday



class PeriodicStats:

    def drop_edges(self, x):
        # should have dont this somewhere up front, bad perfomance of ma is
        # so bad, i want to mask it toward ends
        x[..., 0, :] = np.ma.masked
        x[..., -1, :] = np.ma.masked
        x[..., :, 0] = np.ma.masked
        x[..., :, -1] = np.ma.masked
        return x

    def my_avr(self, x, axis):
        myfnc = [np.average, np.ma.average][x.mask.any().astype(int)]
        o = myfnc(x, axis=axis)
        self.drop_edges(o)
        return o

    def my_med(self, x, axis):
        myfnc = [np.median, np.ma.median][x.mask.any().astype(int)]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            o = myfnc(x, axis=axis)
        self.drop_edges(o)
        return o

    def my_max(self, x, axis):
        myfnc = [np.amax, np.ma.max][x.mask.any().astype(int)]
        o = myfnc(x, axis=axis)
        self.drop_edges(o)
        return o

    def my_cnt(self, x, axis):

        # o = np.ma.MaskedArray(x.count(axis = axis))
        o = np.ma.MaskedArray(x.count(axis=axis), dtype=x.dtype,
                              fill_value=x.fill_value)
        self.drop_edges(o)
        return o

    def __init__(self, fname, oname, period='annual', thres=None):

        self.fname = str(fname)
        if period not in ('annual', 'monthly', 'bimonthly'):
            raise ValueError(f'period: {period}')
        self.period = period

        # get species to process
        spc = [_ for _ in known_spc if _ in self.fname.upper()]
        assert len(spc) == 1
        self.spc = spc[0]
        self.thres = thres

        self.oname = oname

        # pick stats to calculate
        # self.statlist = ['AMAX', 'AAVR', 'AMED', 'A3RD']
        self.statlist = ['AVR', 'MAX', 'MED', 'CNT']

        # description of stats
        # self.statdesc = {k:f'Annual {v}' for k,v in zip(self.statlist,
        #    ['maximum', 'average', 'median', '3rd high'])}
        self.statdesc = {
            'MAX': f'{self.period} maximum',
            'AVR': f'{self.period} mean',
            'MED': f'{self.period} median',
            'CNT': f'{self.period} count of days',
        }
        self.statfunc = {
            #    k:v for k,v in zip(self.statlist,
            # [lambda x: np.amax(x, axis=0),
            #    lambda x: np.average(x, axis=0), 
            #    lambda x: np.median(x, axis=0), 
            #    lambda x: np.percentile(x, q=(363 / 365)*100, axis=0)])}
            'MAX': lambda x: self.my_max(x, axis=0),
            'AVR': lambda x: self.my_avr(x, axis=0),
            'MED': lambda x: self.my_med(x, axis=0),
            'CNT': lambda x: self.my_cnt(x, axis=0),
        }
        self.fo = self.mkheader()
        self.proc()
        self.fo.save(str(self.oname))

    def proc(self):
        for vn0 in self.varlist_in:
            v0 = self.f0.variables[vn0]

            for s in self.statlist:
                vn = s + vn0
                v = self.fo.variables[vn]
                vv = v0[...]
                if self.thres is not None:
                    with np.errstate(invalid='ignore'):
                        vv[v0[...] < self.thres] = np.ma.masked

                for i, (fd0,fd1,tdx0,tdx1) in enumerate(
                        zip(self.firstdays[:-1], self.firstdays[1:], 
                            self.tdx[:-1], self.tdx[1:])):

                    x = self.statfunc[s](vv[tdx0:tdx1, ...])
                    if waste_first_frame:
                        #try:
                        #    v[i + 1, 0, :, :] = x.filled()
                        #except AttributeError:
                        #    v[i + 1, 0, :, :] = x
                        v[i + 1, 0, :, :] = x
                    else:
                        #try:
                        #    v[i, 0, :, :] = x.filled()
                        #except AttributeError:
                        #    v[i, 0, :, :] = x
                        v[i, 0, :, :] = x
                    print(v[0,0,0,0])
                    print(dir(v))


    # make sure that both file follows yosuke's convention to start with
    # last day of previous year
    @staticmethod
    def chktim(f):
        # shoud be daily res
        if f.TSTEP != 240000:
            raise ValueError(f'TSTEP <> 240000: {f.TSTEP}')

        # determine the period covered
        if waste_first_frame:
            # second step should be jan 1 to work around verdi's bug
            dte0 = iodate_to_pydate(f.variables['TFLAG'][1, 0, 0])
            #assert f.variables['TFLAG'][1, 0, 0] % 1000 == 1
        else:
            dte0 = iodate_to_pydate(f.variables['TFLAG'][0, 0, 0])
            #assert f.variables['TFLAG'][0, 0, 0] % 1000 == 1
        dte1 = iodate_to_pydate(f.variables['TFLAG'][-1, 0, 0])

        dates = [dte0 + oned * i
                for i in range((dte1-dte0).days+1)]

        print(dates[-1], dte1)
        assert dates[-1] == dte1
        return dates



    def mkheader(self):

        f0 = pnc.pncopen(self.fname)
        self.f0 = f0
        dates = self.chktim(f0)

        self.year = dates[0].year

        print(dates[0], f0.SDATE // 1000, self.year)
        #if waste_first_frame:
        #    assert f0.SDATE // 1000 - self.year == -1
        #else:
        #    assert f0.SDATE // 1000 - self.year == 0

        self.daysperyear = (datetime.date(self.year + 1, 1, 1) -
                            datetime.date(self.year, 1, 1)).days
        # funky time period to make the timestamp falls approximately the
        # first day of each month
        self.dayspermonth = 30.5
        self.dayspertwomonths = 61

        # list of the begining of each period
        if self.period == 'annual':
            firstday_candidates = [datetime.date(self.year, 1, 1)]
            self.daysperframe = self.daysperyear
        elif self.period == 'monthly':
            firstday_candidates = [datetime.date(self.year, mo, 1) for mo in range(1, 13)]
            self.daysperframe = self.dayspermonth
        elif self.perid == 'bimonthly':
            firstday_candidates = [datetime.date(self.year, mo, 1) for mo in range(2, 13, 2)]
            self.daysperframe = self.dayspertwomonths
        else:
            raise ValueError(f'self.period: {self.period}')
        firstday_candidates.append(datetime.date(self.year+1, 1, 1))

        # firstdays are list of first day of period identified, plus the first day of next period
        # tdx is index for the above
        # so, if for example, one period are found (e.g. asking annual stats for annual dataset), 
        # there would be two elements
        print(firstday_candidates)
        tdx = []
        firstdays = []
        for (fd0, fd1) in zip(firstday_candidates[:-1],firstday_candidates[1:]):
            if fd0 >= dates[0]:
                # first day found
                firstdays.append(fd0)
                if fd1 - oned > dates[-1]:
                    # last day of the period was not found
                    # keep the tdx,firstday, but exit loop
                    tdx.append(len(dates))
                    break
                else:
                    tdx.append(dates.index(fd0))

        if len(tdx) < 2:
            print(firstdays)
            raise ValueError(f'file spans none of specified period: {dates[0]}, {dates[-1]}, ndays={(dates[-1]-dates[0]).days+1}')
        self.firstdays = firstdays
        self.tdx = tdx

        atts = f0.getncatts()

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        statlist = self.statlist
        statdesc = self.statdesc

        # copy the dimensions
        for k, v in f0.dimensions.items():
            if k == 'TSTEP':
                n = len(tdx) - 1
                if waste_first_frame:
                    n += 1
            elif k == 'VAR':
                print(v)
                n = len(statlist)  # max, mean, median, third-high
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # list of output variables
        self.varlist_in = []
        self.varlist_out = []
        for i, nm in enumerate(known_spc):
            if nm not in f0.variables: continue
            self.varlist_in.append(nm)
            for s in statlist:
                self.varlist_out.append(s + nm)

        varlist = self.varlist_out
        # set global attr
        atts['NVARS'] = len(varlist)
        atts['VAR-LIST'] = ''.join([_.ljust(16) for _ in varlist])
        print(f0.dimensions['TSTEP'])
        atts['TSTEP'] = 240000 * self.daysperframe

        sd = self.firstdays[0]
        if waste_first_frame:
            sd -=  oned * self.daysperframe
        atts['SDATE'] = pydate_to_iodate(sd)

        fo.setncatts(atts)

        # make variables
        fo.updatetflag()
        for i, nm in enumerate(['X', 'Y', 'longitude', 'latitude', ]):
            v0 = f0.variables[nm]

            v = fo.createVariable(v0.name, v0.dtype.kind, v0.dimensions)
            v.setncatts(v0.__dict__)
            v[...] = v0[...]

        for nm in self.varlist_in:

            v0 = f0.variables[nm]

            for s in statlist:
                v = fo.createVariable(s + nm, v0.dtype.kind, v0.dimensions)
                a = v0.__dict__
                a['long_name'] = s + a['long_name']
                a['var_desc'] = statdesc[s] + ' ' + a['var_desc']

                v.setncatts(a)

        return fo


def tester():
    Annual_stats(
        'camxv700.MDA8O3.tceq_bc.nc',
        'camxv700.stats.MDA8O3.tceq_bc.nc', )


if __name__ == '__main__':
    tester()
