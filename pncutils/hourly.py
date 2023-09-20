
import PseudoNetCDF as pnc
import numpy as np
import os, re

class Hourly:
    def __init__(self, fnames, oname=None, spc='PM25'):
        """Hourly output joined into one

        :param fnames: list of input file names
        :param oname: (optional) output file name
        :param spc: (optional) output species name
        """

        self.fnames = fnames
        self.spc = spc
        self.raw_spc = [spc]

        not_file = [_ for _ in fnames if not os.path.exists]
        if not_file:
            raise FileNotFoundError('\n'.join(not_file))

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
        i = 0
        for fn in self.fnames:
            # read the next file
            f = pnc.pncopen(fn)
            print(f.SDATE, i)

            for j, s in enumerate(self.raw_spc):
                if j == 0:
                    a = np.array(f.variables[s])
                else:
                    a += np.array(f.variables[s])

            # write
            self.fo.variables[self.spc][i*24:(i+1)*24, ...] = a[...]

            self.fo.updatetflag()

            i += 1

    def _mkheader(self):

        f0 = pnc.pncopen(self.fnames[0])

        atts = f0.getncatts()
        if 'IOAPI_VERSION' not in atts:
            # its not IOAPI file
            raise ValueError(f'not having IOAPI as global attr: {self.fnames[0]}')

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        # copy the dimensions
        for k, v in f0.dimensions.items():
            if k == 'TSTEP':
                n = len(self.fnames) * 24
            elif k == 'VAR':
                n = 1
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # set global attr
        atts['NVARS'] = 1
        atts['VAR-LIST'] = (self.spc).ljust(16)
        fo.setncatts(atts)

        # make variables
        fo.updatetflag()
        for i, nm in enumerate(['X', 'Y', 'longitude', 'latitude', ]):
            v0 = f0.variables[nm]

            v = fo.createVariable(v0.name, v0.dtype.kind, v0.dimensions)
            v.setncatts(v0.__dict__)
            if i < 4:
                v[...] = v0[...]

        v0 = f0.variables[self.raw_spc[0]]
        v = fo.createVariable(self.spc, 'f', v0.dimensions)
        v.setncatts({k: re.sub(self.raw_spc[0], self.spc, v) for k, v in
                     v0.__dict__.items()})
        return fo
