import PseudoNetCDF as pnc
import numpy as np
import re

pm25_species = ['PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4',
                'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'NA', 'PCL']


class Pm25:
    def __init__(self, fname, oname):
        self.fname = fname
        self.oname = oname

        self.fo = self.mkheader()

        self.proc()

        self.fo.save(self.oname)

    def proc(self):
        # read first file

        fn = self.fname

        # self.buf = np.array(pnc.pncopen(fn).variables['O3'])

        i = 0
        # read the next file
        f = pnc.pncopen(self.fname)

        for j, s in enumerate(pm25_species):
            if j == 0:
                a = np.array(f.variables[s])
            else:
                a += np.array(f.variables[s])

        # write
        self.fo.variables['PM25'][...] = a[...]

        self.fo.updatetflag()

    def mkheader(self):

        f0 = pnc.pncopen(self.fname)
        atts = f0.getncatts()

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        # copy the dimensions
        for k, v in f0.dimensions.items():
            if k == 'VAR':
                n = 1
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # set global attr
        atts['NVARS'] = 1
        atts['VAR-LIST'] = 'PM25'.ljust(16)
        fo.setncatts(atts)

        # make variables
        fo.updatetflag()
        for i, nm in enumerate(['X', 'Y', 'longitude', 'latitude', ]):
            v0 = f0.variables[nm]

            v = fo.createVariable(v0.name, v0.dtype.kind, v0.dimensions)
            v.setncatts(v0.__dict__)
            if i < 4:
                v[...] = v0[...]

        v0 = f0.variables['PNH4']
        v = fo.createVariable('PM25', 'f', v0.dimensions)
        v.setncatts({k: re.sub('PNH4', 'A24 PM2.5', v) for k, v in
                     v0.__dict__.items()})
        return fo


if __name__ == '__main__':
    import sys

    fname = sys.argv[1]
    try:
        oname = sys.argv[2]
    except IndexError:
        assert fname[-3:] == '.nc'
        oname = fname[:-3] + '.PM25.nc'
    Pm25(fname, oname)
