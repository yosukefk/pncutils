
import PseudoNetCDF as pnc
import numpy as np
import os, re


class Hourly:

    def from_outputs(self, oname):
        self.oname = oname
        self.fo = pnc.pncopen(self.oname)

    def __init__(self, fnames=None, oname=None, spc='PM25', raw_spc=None):
        """Hourly output joined into one

        :param fnames: list of input file names
        :param oname: (optional) output file name
        :param spc: (optional) output species name
        :param raw_spc: list of species (or function to find species), and sum of them are used as the output species
        """
        self.fnames = fnames
        self.spc = spc
        if not raw_spc:
            self.raw_spc = [spc]
        else:
            self.raw_spc = raw_spc

        if self.fnames is None:
            # bootleg from output
            self._from_outputs(self.oname)
        else:
            self._init(fnames, oname, spc, raw_spc, )


    def _init(self, fnames, oname, spc, raw_spc):
        """Hourly output joined into one

        :param fnames: list of input file names
        :param oname: (optional) output file name
        :param spc: (optional) output species name
        :param raw_spc: list of species (or function to find species), and sum of them are used as the output species
        """

        self.fnames = fnames
        self.spc = spc
        if not raw_spc:
            self.raw_spc = [spc]
        else:
            self.raw_spc = raw_spc
            

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

            
            if callable(self.raw_spc):
                lst = [_ for _ in f.variables if self.raw_spc(_)]
            else:
                lst = self.raw_spc
            print(lst)



            for j, s in enumerate(lst):
                if j == 0:
                    a = np.array(f.variables[s])
                else:
                    a += np.array(f.variables[s])

            # write
            if a.shape[0] < 24: break
            self.fo.variables[self.spc][i*24:(i+1)*24, ...] = a[...]

            self.fo.updatetflag()

            i += 1
        f.close()
        del f

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

        if callable(self.raw_spc):
            lst = [_ for _ in f0.variables if self.raw_spc(_)]
        else:
            lst = self.raw_spc
        v0 = f0.variables[lst[0]]
        v = fo.createVariable(self.spc, 'f', v0.dimensions)
        v.setncatts({k: re.sub(lst[0], self.spc, v) for k, v in
                     v0.__dict__.items()})
        f0.close()
        del f0
        return fo
