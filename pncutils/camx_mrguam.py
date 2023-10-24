import PseudoNetCDF as pnc

class MrgUam:

    def __init__(self, fnames, oname):
        self.fnames = fnames
        self.oname = oname

        self.fo = self._mkheader()
        if self.fo is None:
            return
        print(self.fo.dimensions)


        self._proc()

        if oname:
            self.save(oname)

    def save(self, oname):
        """save

        :param oname: output file name
        """
        self.fo.save(str(oname))

    def _proc(self):
        f = [pnc.pncopen(str(_)) for _ in self.fnames]

        for ifile, ff in enumerate(f):
            for nm, var in ff.variables.items():
                if 'VAR' in var.dimensions:
                    assert var.dimensions[1] == 'VAR'
                    print('\n', nm, var.shape, self.fo.variables[nm].shape, '\n')
                    for j in range(len(self.fo.dimensions['VAR'])):
                        self.fo.variables[nm][:, j, ...] = var[:, 0, ...]

                else:
                    print('\n', nm, var.shape, self.fo.variables[nm].shape, '\n')
                    self.fo.variables[nm][...] = var[...]

    def _mkheader(self):
        # attributes from input files
        f = [pnc.pncopen(str(_)) for _ in self.fnames]
        atts = [_.getncatts() for _ in f]

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        # TODO check compatibility

        # species mapping
        varlists = [
                [att['VAR-LIST'][_*16:(_+1)*16].strip() for _ in range(len(att['VAR-LIST']) // 16)] 
                for att in atts]
        varlists2 = varlists[0].copy()
        for vl in varlists[1:]:
            print(varlists2)
            varlists2 += [_ for _ in vl if _ not in varlists2]
        print(varlists2)

        # create dimensions
        for k,v in f[0].dimensions.items():
            if k == 'VAR':
                n = v.size - len(varlists[0]) + len(varlists2)
            else:
                n = v.size
            fo.createDimension(k, n)
        fo.dimensions['TSTEP'].setunlimited(True)
        #print(fo.dimensions)

        # copy attributes
        atts2 = atts[0].copy()
        atts2['NVARS'] = atts[0]['NVARS'] - len(varlists[0]) + len(varlists2)
        atts2['VAR-LIST'] = ''.join([_.ljust(16) for _ in varlists2])

        fo.setncatts(atts2)


        # make variables
        fo.updatetflag(startdate=f[0].getTimes()[0])

        for ifile, ff in enumerate(f):
            for nm, var in ff.variables.items():
                if nm in fo.variables: continue
                #print(var.name, var.dtype, var.dimensions)

                v = fo.createVariable(var.name, var.dtype, var.dimensions)
                v.setncatts(var.__dict__)


        # createVariable poplulates VAR-LIST, VAR, NVARS
        fo.setncatts({k:v for k,v in atts2.items() if k in ('VAR-LIST', 'NVARS')})

        return fo

def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('filenames',  nargs='+', help='elevated emission file name')
    p.add_argument('-o', '--outname', help='output file name')

    args = p.parse_args()

    ptm = MrgUam(args.filenames, args.outname, )

if __name__ == '__main__':
    Warning('Need QA!')
    main()
