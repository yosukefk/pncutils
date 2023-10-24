# convert elevated emission to lo local format

import PseudoNetCDF as pnc
import warnings


class ElWindow:

    def __init__(self, fname, loname, oname=None, species=None, year=None):
        """camx elevated to lolevel

        :param fname: input elevated file name
        :param loname: reference gridded file name
        :param oname: (optional) output file name
        """

        self.fname = fname
        self.loname = loname
        assert (species is None or isinstance(species, list))
        self.species = species
        self.year = year

        self.fo = self._mkheader()
        if self.fo is None:
            return



        self._proc()

        if oname:
            self.save(oname)

    def save(self, oname):
        """save

        :param oname: output file name
        """
        self.fo.save(oname)

    def _proc(self):
        # attributes from existing files
        f0 = pnc.pncopen(self.fname)
        atts0 = f0.getncatts()
        f1 = pnc.pncopen(self.loname)
        atts1 = f1.getncatts()

        atts = self.fo.getncatts()

        x = f0.variables['xcoord']
        y = f0.variables['ycoord']

        #species = [atts['VAR-LIST'][_*16:(_+1)*16].strip() for _ in range(len(atts['VAR-LIST']) // 16)]
        #for s in species:
        #for s in self.species:
        for s in f0.variables:
            print(s)
            try:
                v0 = f0.variables[s]
            except KeyError:
                continue
            if 'COL' in v0.dimensions:
                icol = v0.dimensions.index('COL')
                print(s)
                print('rhs:', v0)
                print('lhs:',  self.fo.variables[s])
                if icol == 0:
                    self.fo.variables[s][...] = v0[self.myfilter, ...]
                elif icol == 1:
                    self.fo.variables[s][...] = v0[:, self.myfilter, ...]
                elif icol == 2:
                    self.fo.variables[s][...] = v0[:, :, self.myfilter, ...]
                else:
                    raise
            else:
                self.fo.variables[s][...] = v0[...]





        

    def _mkheader(self):

        # attributes from existing files
        f0 = pnc.pncopen(str(self.fname))
        atts0 = f0.getncatts()
        f1 = pnc.pncopen(str(self.loname))
        atts1 = f1.getncatts()


        # find how many stack are in new domain
        x = f0.variables['xcoord']
        y = f0.variables['ycoord']

        myfilter = (x[...] > f1.XORIG) & (x[...] > (f1.XORIG+f1.XCELL*f1.NCOLS)) & (y[...]  > f1.YORIG) & (y[...]  > (f1.YORIG+f1.YCELL*f1.NROWS))
        print(f'nstk: orig={len(myfilter)}, filtered={myfilter.sum()}')
        nstk = myfilter.sum()

        self.myfilter = myfilter
        self.nstk = nstk

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        # copy the dimensions
        for k, v in f0.dimensions.items():
            if k == 'COL':
                n = nstk
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # copy attributes
        atts = atts0.copy()
        atts.update({k:v for k,v in atts1.items() if k in 
            ('XORIG', 'YORIG', 'XCELL', 'YCELL', 'NCOLS', 'NROWS',
                    )})

        if self.year is not None:
            print(atts['SDATE'])
            atts['SDATE'] = atts['SDATE'] % 1000 + self.year * 1000
            print(atts['SDATE'])
        fo.setncatts(atts)


        # make variables
        fo.updatetflag()

        myspecies = [atts['VAR-LIST'][_*16:(_+1)*16].strip() for _ in range(len(atts['VAR-LIST']) // 16)]
        if self.species is None:
            self.species = myspecies
        else:
            found = [_ for _ in self.species if _ in myspecies]
            if len(found) == 0:
                warnings.warn(f'none of {self.species} found in file')
                return None
            self.species = found
            
        names = self.species
        print(names)

        vx = f1.variables[atts1['VAR-LIST'][:16].strip()]

        #for i,nm in enumerate(names):
        for i,nm in enumerate(f0.variables):
            if nm == 'TFLAG':
                fo.updatetflag()
                continue

            vv = f0.variables[nm]
            dim = vx.dimensions

            v = fo.createVariable(vv.name, vv.dtype.kind, vv.dimensions)#, compression='zlib', complevel=2)
            v.setncatts(vv.__dict__)

        attso = fo.getncatts()
        return fo




def tester():

    fin = 'camx_cb6p_ei_el_neguus.2019_day.txo3.bc19.2019_pig_v2.nc'
    flo = 'camx_cb6_ei_lo_loptus.2019_day.txo3.bc19.negu_2019_v2.txs_4km.nc'
    fout = 'elwindow.nc'

    lo = ElWindow(fin, flo, fout)
    return lo

#lo = tester()

def main():
    import sys
    import argparse
    import re
    p = argparse.ArgumentParser()
    p.add_argument('filename',  help='elevated emission file name')
    p.add_argument('outname', nargs='?', help='output file name')

    p.add_argument('-s', '--species', help='species to extract (comma or space delimited)')
    p.add_argument('-l', '--loname', required=True, help='lolevel file name (to define grid)')

    p.add_argument('-y', '--year', help='overwrite year', type=int )
    args = p.parse_args()

    if args.species is not None:
        args.species = re.split('[, ] *', args.species)

    if args.outname is None:
        args.outname = args.filename[:-3] + '.elwindow.nc'


    lo = ElWindow(args.filename, args.loname, args.outname, args.species, year=args.year)


if __name__ == '__main__':
    Warning('Need QA!')
    main()



