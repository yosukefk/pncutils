#!/usr/bin/env python
import PseudoNetCDF as pnc
import numpy as np

class Window:


    def __init__(self, fname, loname, oname=None, species=None):
        self.fname = fname
        self.oname = oname
        self.loname = loname

        assert (species is None or isinstance(species, list))
        self.species = species

        self.fo = self._mkheader()
        if self.fo is None:
            return
        #print(self.fo.dimensions)


        self._proc()

        if oname:
            self.save(oname)

    def save(self, oname):
        """save

        :param oname: output file name
        """
        self.fo.save(str(oname))

    def _proc(self):
        f0 = pnc.pncopen(str(self.fname))
        atts0 = f0.getncatts()
        f1 = pnc.pncopen(self.loname)
        atts1 = f1.getncatts()

        atts = self.fo.getncatts()

        for s in self.species:
            print(s)
            try:
                v0 = f0.variables[s]
            except KeyError:
                continue

            if self.meshfac > 1:
                raise NotImplementedError(f'meshfac > 1: {self.meshfac}') 
            else:
                val = v0[:, 0, self.joff:(self.joff+f1.NROWS//self.rat), 
                        self.ioff:(self.ioff+f1.NCOLS//self.rat)]
                val2 = np.repeat(
                        np.repeat(val, self.rat, axis=-1),
                        self.rat, axis=-2)
                val2 *= (self.meshfac * self.meshfac)


            self.fo.variables[s][:, 0, :, :] = val2

    def _mkheader(self):
        # attributes from input files
        f0 = pnc.pncopen(str(self.fname))
        atts0 = f0.getncatts()

        f1 = pnc.pncopen(str(self.loname))
        atts1 = f1.getncatts()

        # mesh factor
        xc0 = f0.XCELL
        xc1 = f1.XCELL

        if xc0 >= xc1:
            rat = xc0 / xc1
        else:
            rat = xc1 / xc0

        irat = round(rat)
        if abs((rat - irat ) / rat) > 0.01:
            raise RuntimeError('uneven ratio:', xc0, xc1, rat)
        meshfac = xc1 / xc0

        # offset
        xo0, yo0 = f0.XORIG, f0.YORIG
        xo1, yo1 = f1.XORIG, f1.YORIG

        xoff = xo1 - xo0
        yoff = yo1 - yo0

        rioff = xoff / xc0
        ioff = round(rioff)
        if abs((rioff - ioff) / rioff) > 0.01:
            raise RuntimeError('uneven xoff:', xo0, xo1, rioff)

        rjoff = yoff / xc0
        joff = round(rjoff)
        if abs((rjoff - joff) / rjoff) > 0.01:
            raise RuntimeError('uneven yoff:', yo0, yo1, rjoff)

        self.meshfac = meshfac
        self.rat = irat
        self.ioff = ioff
        self.joff = joff

        print(f'meshfac: {meshfac}\nrat: {rat}\nioff: {ioff}\njoff: {joff}')


        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        # copy the dimensions
        for k, v in f0.dimensions.items():
            if k in ('COL', 'ROW'):
                n = f1.dimensions[k].size
            elif k == 'VAR':
                n = atts0['NVARS'] + 6
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # copy attributes
        atts = atts0.copy()
        atts.update({k:v for k,v in atts1.items() if k in 
            ('NCOLS', 'NROWS', 'XORIG', 'YORIG', 'XCELL', 'YCELL', 
                    'CAMx_NAME', 'FILEDESC')})
        fo.setncatts(atts)


        # make variables
        fo.updatetflag()
        #names = ['X', 'Y', 'TFLAG', 'ETFLAG', 'longitude', 'latitude'] 
        #names = ['X', 'Y',                    'longitude', 'latitude'] 
        names0 = ['X', 'Y', 'TFLAG', 'ETFLAG', 'longitude', 'latitude'] 

        myspecies = [atts['VAR-LIST'][_*16:(_+1)*16].strip() for _ in range(len(atts['VAR-LIST']) // 16)]
        if self.species is None:
            self.species = myspecies
        else:
            found = [_ for _ in self.species if _ in myspecies]
            if len(found) == 0:
                warnings.warn(f'none of {self.species} found in file')
                return None
            self.species = found
            
        names = names0 + self.species
        print(names)

        vx = f1.variables[atts1['VAR-LIST'][:16].strip()]

        for i,nm in enumerate(names):
            if nm == 'TFLAG':
                fo.updatetflag()
                continue

            if nm in names0:
                vv = f1.variables[nm]
                dim = vv.dimensions
            else:
                vv = f0.variables[nm]
                dim = vx.dimensions

            v = fo.createVariable(vv.name, vv.dtype.kind, dim)#, compression='zlib', complevel=2)
            v.setncatts(vv.__dict__)

            if nm in names0:
                print(nm)
                if nm=='ETFLAG':
                    print(v.shape, vv.shape)
                    v[...] = vv[:, 0, :].repeat(v.shape[1], axis=1).reshape(v.shape)
                    #v[...] = vv[:, :v.shape[1], :]
                else:
                    v[...] = vv[...]
        attso = fo.getncatts()
        return fo

def main():
    import sys
    import argparse
    import re
    p = argparse.ArgumentParser()
    p.add_argument('filename',  help='gridded emission file name')
    p.add_argument('outname', nargs='?', help='output file name')

    p.add_argument('-s', '--species', help='species to extract (comma or space delimited)')
    p.add_argument('-l', '--loname', required=True, help='lolevel file name (to define grid)')

    args = p.parse_args()

    if args.species is not None:
        args.species = re.split('[, ] *', args.species)

    if args.outname is None:
        args.outname = args.filename[:-3] + '.window.nc'

    wn = Window(args.filename, args.loname, args.outname)

if __name__ == '__main__':
    main()
