#!/usr/bin/env python3

import PseudoNetCDF as pnc
import numpy as np
import datetime
import warnings

known_spc = ['A24PM25', 'MDA8O3']

# fudge to work around bug for verdi...
waste_first_frame = True


class diff_stats:

    
    def drop_edges(self, x):
        # should have done this somewhere up front, bad perfomance of ma is
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

    def my_min(self, x, axis):
        myfnc = [np.amin, np.ma.min][x.mask.any().astype(int)]
        o = myfnc(x, axis=axis)
        self.drop_edges(o)
        return o

    def my_max(self, x, axis):
        myfnc = [np.amax, np.ma.max][x.mask.any().astype(int)]
        o = myfnc(x, axis=axis)
        self.drop_edges(o)
        return o

    def my_mxdec(self, x, axis):
        myfnc = [np.amin, np.ma.min][x.mask.any().astype(int)]
        o = myfnc(x, axis=axis)
        o = np.minimum(o, 0)
        self.drop_edges(o)
        return o

    def my_mxinc(self, x, axis):
        myfnc = [np.amax, np.ma.max][x.mask.any().astype(int)]
        o = myfnc(x, axis=axis)
        o = np.maximum(o, 0)
        self.drop_edges(o)
        return o

    def my_cnt(self, x, axis):

        #o = np.ma.MaskedArray(x.count(axis = axis))
        o = np.ma.MaskedArray(x.count(axis = axis), dtype=x.dtype,
                fill_value = x.fill_value)
        self.drop_edges(o)
        return o

    def masked_minmax(self, x, axis):
        if not x.mask.any():
            return self.minmax(x, axis)
        # report min of max, whichever has greater abs value
        a = x.squeeze()
        mx = np.expand_dims(np.ma.max(a, axis=axis), axis=axis)
        mn = np.expand_dims(np.ma.min(a, axis=axis), axis=axis)
        mnmx = np.ma.concatenate((mn,mx), axis=axis)
        # mn, mx has np.nan, and causes wrnings in expression below, so...
        with np.errstate(invalid='ignore'):
            whch = (np.abs(mn) < np.abs(mx)).astype(int)
        o = np.take_along_axis(
                mnmx,
                whch,
                axis=axis
                )
        self.drop_edges(o)
        return o

    def minmax(self, x, axis):
        # report min of max, whichever has greater abs value
        a = x.squeeze()
        mx = np.expand_dims(np.amax(a, axis=axis), axis=axis)
        mn = np.expand_dims(np.amin(a, axis=axis), axis=axis)
        mnmx = np.concatenate((mn,mx), axis=axis)
        with np.errstate(invalid='ignore'):
            whch = (np.abs(mn) < np.abs(mx)).astype(int)
        o = np.take_along_axis(
                mnmx,
                whch,
                axis=axis
                )
        self.drop_edges(o)
        return o

    def __init__(self, fname1, fname0, oname, period = 'annual', thres=None):
        self.fname0 = fname0
        self.fname1 = fname1
        if period not in ('annual', 'monthly', 'bimonthly'):
            raise ValueError(f'period: {period}')
        self.period = period

        # get species to process
        spc = [_ for _ in known_spc if _ in fname1]
        assert len(spc) == 1
        self.spc = spc[0]
        self.thres = thres

        self.oname = oname

        # pick stats to calculate
#        self.statlist = ['MMX','AVR',]
        self.statlist = ['MMXD','AVRD','MAXD','MIND','MEDD','CNTD','MXINC','MXDEC']

        # description of stats
        self.statdesc = {
                'MMXD': f'{self.period} absolute maximum diff', 
                'MAXD': f'{self.period} (positive) maximum diff',
                'MIND': f'{self.period} minimum ("negative max") diff',
                'AVRD': f'{self.period} mean diff',
                'MEDD': f'{self.period} median diff',
                'CNTD': f'{self.period} count of days',
                'MXINC': f'{self.period} maximum increase',
                'MXDEC': f'{self.period} maximum decrease',
                }
        self.statfunc = { 
            # at first it was like this, simple
            #'MAXD': lambda x: np.amax(x, axis=0),
            #'MIND': lambda x: np.amin(x, axis=0),
            # switches depending on any cell got masked
            #'AVRD': lambda x: [np.average, np.ma.average][x.mask.any().astype(int)](x, axis=0),
            #'MEDD': lambda x: [np.median, np.ma.median][x.mask.any().astype(int)](x, axis=0),
            #'CNTD': lambda x: x.count(axis=0),
            # i figured better make these function...  
            'MMXD': lambda x: self.masked_minmax(x, axis=0),
            'MAXD': lambda x: self.my_max(x, axis=0),
            'MIND': lambda x: self.my_min(x, axis=0),
            'AVRD': lambda x: self.my_avr(x, axis=0),
            'MEDD': lambda x: self.my_med(x, axis=0),
            'CNTD': lambda x: self.my_cnt(x, axis=0),
            'MXINC': lambda x: self.my_mxinc(x, axis=0),
            'MXDEC': lambda x: self.my_mxdec(x, axis=0),
                } 
        self.fo = self.mkheader()
        self.proc()
        self.fo.save(self.oname)

    def proc(self):
        # debugging...
        self.o = {}
        self.oo = {}
        for vn0 in self.varlist_in:
            v0 = self.f0.variables[vn0]
            v1 = self.f1.variables[vn0]
            print('v0,v1')
            print( np.ma.max(v0[...]))
            print( np.ma.max(v1[...]))
            for s in self.statlist:
                vn=s+vn0
                v = self.fo.variables[vn]
                vdiff= v1[...] - v0[...]
                assert isinstance(vdiff, np.ma.MaskedArray)
                print('vdiff')
                print(np.ma.max(vdiff))
                if self.thres is not None:
                    #vdiff = np.ma.masked_where(v0[...]< self.thres, vdiff)
                    with np.errstate(invalid='ignore'):
                        vdiff[v0[...] < self.thres] = np.ma.masked
                    
                if self.period == 'annual':
                    if waste_first_frame: 
                        #v[1,0,:,:] = self.statfunc[s](vdiff[1:,...])
                        x = self.statfunc[s](vdiff[1:,...])
                        #self.o[vn] = x
                        try:
                            v[1,0,:,:] = x.filled()
                        except AttributeError:
                            v[1,0,:,:] = x
                        #self.oo[vn] = v[...]
                    else:
                        x = self.statfunc[s](vdiff)
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
                            x = self.statfunc[s](vdiff[day1[i]:day1[i+1],...])
                            try:
                                v[i+1,0,:,:] = x.filled()
                            except AttributeError:
                                v[i+1,0,:,:] = x
                        else:
                            x = self.statfunc[s](vdiff[(day1[i]-1):(day1[i+1]-1),...])
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
                            x = self.statfunc[s](vdiff[day1[i]:day1[i+1],...])
                            try:
                                v[i+1,0,:,:] = x.filled()
                            except AttributeError:
                                v[i+1,0,:,:] = x

                        else:
                            x = self.statfunc[s](vdiff[(day1[i]-1):(day1[i+1]-1),...])
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
        f0 = pnc.pncopen(self.fname0)
        f1 = pnc.pncopen(self.fname1)
        self.f0 = f0
        self.f1 = f1



        self.chktim(f0)
        self.chktim(f1)
        if waste_first_frame:
            self.year = int(f1.variables['TFLAG'][1,0,0] / 1000)
            assert int(f1.SDATE  / 1000) - self.year == -1
        else:
            self.year = int(f1.variables['TFLAG'][0,0,0] / 1000)
            assert int(f1.SDATE  / 1000) - self.year == 0
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

        atts = f1.getncatts()

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
            #atts['SDATE'] = f0.variables['TFLAG'][1,0,0]
            #atts['STIME'] = f0.variables['TFLAG'][1,0,1]
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
    if False:
        o = diff_stats( 
                'camxv700.MDA8O3.v2b_no_egu_coal.nc', 
                'camxv700.MDA8O3.crisp_bc_v2b.nc', 
                'test_camxv700.annualdiff.MDA8O3.v2b_no_egu_coal_minus_crisp_bc_v2.nc'
                )
    if False:
        o = diff_stats( 
                'camxv700.MDA8O3.v2b_no_egu_coal.nc', 
                'camxv700.MDA8O3.crisp_bc_v2b.nc', 
                'test_camxv700.monthlydiff.MDA8O3.v2b_no_egu_coal_minus_crisp_bc_v2.nc', 
                period='monthly',
                )
    if True:
        o = diff_stats( 
                'camxv700.MDA8O3.v2b_no_egu.nc', 
                'camxv700.MDA8O3.crisp_bc_v2b.nc', 
                'test_camxv700.annualdiff.MDA8O3_thres60.v2b_no_egu_minus_crisp_bc_v2.nc', 
                thres = 0.06
                )
        return o

if __name__ == '__main__':
    tester(

        )

