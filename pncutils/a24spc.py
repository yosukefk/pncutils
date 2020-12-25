import PseudoNetCDF as pnc

import numpy as np
import re

pm25 = ['PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4',
        'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'NA', 'PCL']


class A24spc:
    def __init__(self, fnames, oname=None):
        """24 hour average PM2.5, unaggregated

        :param fnames: list of input file names
        :param oname: (optional) output file name
        """
        self.fnames = fnames

        self.fo = self._mkheader()

        self._proc()

        if oname:
            self.save(oname)

    def save(self, oname):
        """save

        :param oname: output file name
        """
        self.fo.save(oname)

    def _proc(self):
        # read first file

        fn = self.fnames[0]

        # self.buf = np.array(pnc.pncopen(fn).variables['O3'])

        i = 0
        for fn in self.fnames:
            # read the next file
            f = pnc.pncopen(fn)
            print(f.SDATE, i)

            for j, s in enumerate(pm25):
                a = np.array(f.variables[s])

                a24 = a.mean(axis=0)

                # get rid of values at boudary cells
                a24[..., 0, :] = np.nan
                a24[..., -1, :] = np.nan
                a24[..., :, 0] = np.nan
                a24[..., :, -1] = np.nan

                # write
                self.fo.variables[s][i, ...] = a.mean(axis=0)

            self.fo.updatetflag()

            i += 1

    def _mkheader(self):

        f0 = pnc.pncopen(self.fnames[0])
        atts = f0.getncatts()

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        # copy the dimensions
        for k, v in f0.dimensions.items():
            if k == 'TSTEP':
                n = len(self.fnames)
            elif k == 'VAR':
                n = 1
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # set global attr
        atts['NVARS'] = len(pm25)
        atts['VAR-LIST'] = ''.join(_.ljust(16) for _ in pm25)
        atts['TSTEP'] = 240000
        fo.setncatts(atts)

        # make variables
        fo.updatetflag()
        for i, nm in enumerate(['X', 'Y', 'longitude', 'latitude', ]):
            v0 = f0.variables[nm]

            v = fo.createVariable(v0.name, v0.dtype.kind, v0.dimensions)
            v.setncatts(v0.__dict__)
            if i < 4:
                v[...] = v0[...]

        for s in pm25:
            v0 = f0.variables[s]
            v = fo.createVariable(s, 'f', v0.dimensions)
            v.setncatts({k: v for k, v in v0.__dict__.items()})
        return fo


def tester():
    o3 = A24spc(
        """
camx65ss_cb6r4hCF.20160101.rh.bc16_16s1.v1LN_v1a.2016_wrf381_p2KFsn_i2KFsn.avrg.grd01.nc
camx65ss_cb6r4hCF.20160102.rh.bc16_16s1.v1LN_v1a.2016_wrf381_p2KFsn_i2KFsn.avrg.grd01.nc
camx65ss_cb6r4hCF.20160103.rh.bc16_16s1.v1LN_v1a.2016_wrf381_p2KFsn_i2KFsn.avrg.grd01.nc
""".strip().split()
        , 'ooo.nc')
    return o3

# o3 = tester()
