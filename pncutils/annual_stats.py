#!/usr/bin/env python3

import PseudoNetCDF as pnc
import numpy as np
import datetime
import warnings

known_spc = ['A24PM25', 'MDA8O3']

# fudge to work around bug for verdi...
waste_first_frame = True


class annual_stats:


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

        #o = np.ma.MaskedArray(x.count(axis = axis))
        o = np.ma.MaskedArray(x.count(axis = axis), dtype=x.dtype,
                fill_value = x.fill_value)
        self.drop_edges(o)
        return o

    def __init__(self, fname, oname, period = 'annual', thres=None):

        self.fname = fname
        if period not in ('annual', 'monthly', 'bimonthly'):
            raise ValueError(f'period: {period}')
        self.period = period

        # get species to process
        spc = [_ for _ in known_spc if _ in fname]
        assert len(spc) == 1
        self.spc = spc[0]
        self.thres = thres

        self.oname = oname

        # pick stats to calculate
        #self.statlist = ['AMAX', 'AAVR', 'AMED', 'A3RD']
        self.statlist = ['AVR','MAX','MED','CNT']

        # description of stats
        #self.statdesc = {k:f'Annual {v}' for k,v in zip(self.statlist,
        #    ['maximum', 'average', 'median', '3rd high'])}
        self.statdesc = {
                'MAX': f'{self.period} maximum',
                'AVR': f'{self.period} mean',
                'MED': f'{self.period} median',
                'CNT': f'{self.period} count of days',
                }
        self.statfunc = {
            #    k:v for k,v in zip(self.statlist,
            #[lambda x: np.amax(x, axis=0), 
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
        self.fo.save(self.oname)

    def proc(self):
        for vn0 in self.varlist_in:
            v0 = self.f0.variables[vn0]
            for s in self.statlist:
                vn=s+vn0
                v = self.fo.variables[vn]
                vv = v0[...]
                if self.thres is not None:
                    with np.errstate(invalid='ignore'):
                        vv[v0[...] < self.thres] = np.ma.masked

                if self.period == 'annual':
                    if waste_first_frame: 
                        x = self.statfunc[s](vv[1:,...])
                        try:
                            v[1,0,:,:] = x.filled()
                        except AttributeError:
                            v[1,0,:,:] = x
                    else:
                        x = self.statfunc[s](vv)
                        try:
                            v[0,0,:,:] = x.filled()
                        except AttributeError:
                            v[0,0,:,:] = x

                elif self.period == 'monthly':

                    # first day of each month
                    day1 = [datetime.date(self.year,(_+1),1) for _ in
                            range(12)]
                    # append first day of next year
                    day1.append(datetime.date(self.year+1,1,1))

                    # express the date as jday
                    day1 = [(_ - day1[0]).days + 1 for _ in day1]

                    for i in range(12):
                        if waste_first_frame: 
                            # since 12/31 of previous year is included, day1 can be used
                            # directly as index of python array
                            x = self.statfunc[s](vv[day1[i]:day1[i+1],...])
                            try:
                                v[i+1,0,:,:] = x.filled()
                            except AttributeError:
                                v[i+1,0,:,:] = x
                        else:
                            x = self.statfunc[s](vv[(day1[i]-1):(day1[i+1]-1),...])
                            try:
                                v[i,0,:,:] = x.filled()
                            except AttributeError:
                                v[i,0,:,:] = x
                elif self.period == 'bimonthly':
                    # first day of odd month
                    day1 = [datetime.date(self.year,(2*_+1),1) for _ in
                            range(6)]
                    # append first day of next year
                    day1.append(datetime.date(self.year+1,1,1))
                    # express the date as jday
                    day1 = [(_ - day1[0]).days + 1 for _ in day1]
                    for i in range(6):
                        if waste_first_frame: 
                            # since there is 12/31 included, day1 can be used
                            # directly as index of python array
                            x = self.statfunc[s](vv[day1[i]:day1[i+1],...])
                            try:
                                v[i+1,0,:,:] = x.filled()
                            except AttributeError:
                                v[i+1,0,:,:] = x

                        else:
                            x = self.statfunc[s](vv[(day1[i]-1):(day1[i+1]-1),...])
                            try:
                                v[i,0,:,:] = x.filled()
                            except AttributeError:
                                v[i,0,:,:] = x
                else:
                    raise ValueError(f'self.period: {self.period}')

    # make sure that both file follows yosuke's convention to start with
    # last day of previous year
    @staticmethod
    def chktim(f):
        # shoud be daily res
        assert f.TSTEP == 240000
        if waste_first_frame:
            # second step should be jan 1 to work around verdi's bug
            assert f.variables['TFLAG'][1,0,0] % 1000 == 1
        else:
            assert f.variables['TFLAG'][0,0,0] % 1000 == 1
        # last step shouls be last day
        assert (f.variables['TFLAG'][-1,0,0] % 1000) in (365,366)


    def mkheader(self):
        f0 = pnc.pncopen(self.fname)
        self.f0 = f0
        self.chktim(f0)
        if waste_first_frame:
            self.year = int(f0.variables['TFLAG'][1,0,0] / 1000)
            assert int(f0.SDATE  / 1000) - self.year == -1
        else:
            self.year = int(f0.variables['TFLAG'][0,0,0] / 1000)
            assert int(f0.SDATE  / 1000) - self.year == 0
        self.daysperyear = (datetime.date(self.year+1, 1, 1) -
                datetime.date(self.year, 1, 1)).days
        # funky time period to make the timestamp falls approximately the
        # first day of each month
        self.dayspermonth = 30.5 
        self.dayspertwomonths = 61

        if self.period == 'annual':
            self.daysperframe = self.daysperyear
        elif self.period == 'monthly':
            self.daysperframe = self.dayspermonth
        elif self.period == 'bimonthly':
            self.daysperframe = self.dayspertwomonths
        else:
            raise ValueError(f'self.period: {self.period}')

        atts = f0.getncatts()

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        statlist = self.statlist
        statdesc = self.statdesc

        # copy the dimensions
        for k,v in f0.dimensions.items():
            if k == 'TSTEP':
                if self.period == 'annual':
                    n = 1
                elif self.period == 'monthly':
                    n = 12
                elif self.period == 'bimonthly':
                    n = 6
                else:
                    raise ValueError(f'self.period: {self.period}')

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
                self.varlist_out.append(s+nm)


        varlist = self.varlist_out
        # set global attr
        atts['NVARS'] = len(varlist)
        atts['VAR-LIST'] = ''.join([_.ljust(16) for _ in varlist])
        print(f0.dimensions['TSTEP'])
        atts['TSTEP'] = 240000 * self.daysperframe

        if waste_first_frame:
            sd = datetime.datetime(self.year, 1, 1) - datetime.timedelta(
                    days=self.daysperframe)
            atts['SDATE'] = sd.year * 1000 + ( (sd - datetime.datetime(sd.year, 1,
                1)).days + 1 )
            atts['STIME'] = sd.hour * 10000 + sd.minute * 100 + sd.second
        else:
            pass
        fo.setncatts(atts)

        # make variables
        fo.updatetflag()
        for i,nm in enumerate(['X', 'Y', 'longitude', 'latitude', ]):
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
                a['var_desc'] = statdesc[s] + ' ' +a['var_desc']
                
                v.setncatts(a)



        return fo

def tester():
    annual_stats(
'camxv700.MDA8O3.tceq_bc.nc',
'camxv700.stats.MDA8O3.tceq_bc.nc',)

if __name__ == '__main__':
    tester()

