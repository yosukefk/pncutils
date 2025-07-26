import hourly as hrly
import PseudoNetCDF as pnc
import numpy as np
from pathlib import Path
import re


class Daily:

    def __init__(self, fnames=None, hourly_name=None, oname=None, 
            spc='MDA8O3', hourly_spc=None, raw_spc=None, fnc_dayagg=None, itzon_use=None):

        # grab data
        if hourly_name is None:
            if fnames is None:
                raise RuntimeError('neither hourly_name or fnames is spcified')
            self.fnames = fnames
        else:
            if Path(hourly_name).is_file():
                # read from preprocessed hourly file
                self.hourly = pnc.pncopen(hourly_name)
            else:
                # need to generate hourly dataset
                if fnames is None:
                    raise RuntimeError('hourly_name doesnt exist and fnames is spcified')
                else:
                    self.hourly = hrly.Hourly(fnames, oname=None, spc=spc, raw_spc=raw_spc)

        # timezone
        if itzon_use is None: itzon_use = self.hourly.ITZON
        self.itzon_use = itzon_use
        print(self.hourly.ITZON, itzon_use)
        self.tzoffset = self.hourly.ITZON - itzon_use

        self.spc = spc
        self.hourly_spc = hourly_spc
        self.fnc_dayagg = fnc_dayagg

        self.fo = self._mkheader()
        self._proc()
        if oname:
            self.save(oname)

    def _mkheader(self):
        atts = self.hourly.getncatts()
        if 'IOAPI_VERSION' not in atts:
            # its not IOAPI file
            raise ValueError(f'not having IOAPI as global attr: {self.hourly}')

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()


        # copy the dimensions
        for k, v in self.hourly.dimensions.items():
            if k == 'TSTEP':
                self.ndays = (v.size - abs(self.tzoffset)) // 24
                n = self.ndays
            elif k == 'VAR':
                n = 1
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # set global attr
        atts['NVARS'] = 1
        atts['VAR-LIST'] = (self.spc).ljust(16)
        atts['TSTEP'] = 240000
        atts['ITZON'] = self.itzon_use
        fo.setncatts(atts)

        # make variables
        fo.updatetflag()
        for i, nm in enumerate(['X', 'Y', 'longitude', 'latitude', ]):
            v0 = self.hourly.variables[nm]

            v = fo.createVariable(v0.name, v0.dtype.kind, v0.dimensions)
            v.setncatts(v0.__dict__)
            if i < 4:
                v[...] = v0[...]

        v0 = self.hourly.variables[self.hourly_spc]
        v = fo.createVariable(self.spc, 'f', v0.dimensions)
        v.setncatts({k: re.sub(self.hourly_spc, self.spc, v) for k, v in
                     v0.__dict__.items()})
        return fo

    def _proc(self):
        if self.tzoffset <= 0:
            self.buf = self.hourly.variables[self.hourly_spc][-self.tzoffset:, ...]
            dayoffset = 0
        else:
            self.buf = self.hourly.variables[self.hourly_spc][24-self.tzoffset:, ...]
            self.fo.variables[self.spc][0, ...] = np.nan
            self.fo.updatetflag()
            dayoffset = 1
            raise RuntimeError('need to write code, having np.nam as first tstep')
        out = self.fnc_dayagg(self.buf)

        for i in range(self.ndays):
            self.fo.variables[self.spc][(i + dayoffset), ...] = out[i, ...]
            self.fo.updatetflag()



    def save(self,oname):
        """save

        :param oname: output file name
        """
        self.fo.save(oname)


class MDA8O3(Daily):
    def __init__(self, fnames=None, hourly_name=None, oname=None, itzon_use=None):
        super().__init__(fnames=fnames, hourly_name=hourly_name, oname=oname,  
                itzon_use=itzon_use, spc='MDA8O3', hourly_spc='O3', fnc_dayagg=calc_mda8)

class A24PM25(Daily):
    def __init__(self, fnames=None, hourly_name=None, oname=None, itzon_use=None):
        super().__init__(fnames=fnames, hourly_name=hourly_name, oname=oname,  
                itzon_use=itzon_use, spc='A24PM25', hourly_spc='PM25', fnc_dayagg=calc_a24) 
        
def calc_a24(x):
    shp = list(x.shape)
    shp[0] = ( shp[0]  - 1 ) // 24 + 1
    out = np.empty(shp) 
    out[-1, ...] = np.nan
    for i in range(shp[0]):
        out[i, ...] = a8[(24*i):(24*(i+1))].mean(axis=0)
    return out


def calc_mda8(x):
    a8 = moving_average(x, 8)
    shp = list(x.shape)
    shp[0] = ( shp[0]  - 1 ) // 24 + 1
    out = np.empty(shp) 
    out[-1, ...] = np.nan
    for i in range(shp[0]):
        out[i, ...] = a8[(24*i):(24*(i+1))].max(axis=0)
    return out


def moving_average(x, w=8):
    """moving average on the first axis

    :param x: numpy float array
    :param w: (optional) window size, default 8
    :return: numpy float array, first dim shorter py w-1
    """
    n = x.shape[0]
    m = n - w + 1
    o = np.zeros((m,) + x.shape[1:])
    for i in range(n):
        for j in range(w):
            if i - j < 0 or i - j >= m: continue
            o[(i - j), ...] += x[i, ...]
    return o / w
