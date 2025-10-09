# merge gridded emission files


import PseudoNetCDF as pnc
import pandas as pd, numpy as np
import warnings


class MrgUam:
    def __init__(self, fnames, oname, summary=None):
        self.fnames = fnames
        self.oname = oname
        self.summary = []

        self.fo = self._mkheader()
        if self.fo is None:
            return
        #print(self.fo.dimensions)


        self._proc()

        if oname:
            self.save(oname)
        if summary is not None:
            df_summary=pd.DataFrame(self.summary)
            df_summary.to_csv(summary)
            

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
                    #print('\n', nm, var.shape, self.fo.variables[nm].shape, '\n')

                    if ifile == 0 or nm in ('TFLAG', 'ETFLAG'):
                        # if first file, just copy
                        # or for TFLAG, simply overwrite with each file's value (last file wins)
                        for j in range(len(self.fo.dimensions['VAR'])):
                            self.fo.variables[nm][:, j, ...] = var[:, 0, ...]
                    else:
                        # check if values are consistent
                        for j in range(len(self.fo.dimensions['VAR'])):
                            try:
                                assert np.array_equal(self.fo.variables[nm][:, j, ...] , var[:, 0, ...])
                            except AssertionError:
                                print(nm, j,)
                                raise

                else:
                    #print('\n', nm, var.shape, self.fo.variables[nm].shape, '\n')
                    if nm in self.varlists2:
                        self.fo.variables[nm][...] += var[...]

                        v = self.fo.variables[nm]
                        if v.units.endswith('hr-1') and v.shape[0] == 25:
                            mysum = var[:24, ...].sum()
                        else:
                            mysum = var[...].sum()

                        self.summary.append({
                            'fname': self.fnames[ifile],
                            'spc': nm,
                            'shp': var.shape,
                            'units': var.units,
                            'sum': mysum,
                            })
                    else:
                        if ifile == 0:
                            warnings.warn(f'not in varlists2: {nm}, fname = {self.fnames[ifile]}')
                            self.fo.variables[nm][...] = var[...]

                            v = self.fo.variables[nm]
                            if v.units.endswith('hr-1') and v.shape[0] == 25:
                                mysum = var[:24, ...].sum()
                            else:
                                mysum = var[...].sum()

                            self.summary.append({
                                'fname': self.fnames[ifile],
                                'spc': nm,
                                'shp': var.shape,
                                'units': var.units,
                                'sum': mysum,
                                })
                        else:
                            assert np.array_equal(self.fo.variables[nm][...] , var[...])

    def _mkheader(self):
        # attributes from input files
        f = [pnc.pncopen(str(_)) for _ in self.fnames]
        atts = [_.getncatts() for _ in f]

        # check compatibility of TSTEP
        ntims = [len(_.dimensions['TSTEP']) for _ in f]
        self.ignore_last_tstep = False
        if all(_ == ntims[0] for _ in  ntims):
            # ok
            pass
        elif all(_ in (24, 25) for _ in ntims) :
            # special case
            # only if files are either 24 and 25, and time of day are 0-23 or 0-24
            #tflags = [np.array(_.variables['TFLAG'][:,0,1]) for _ in f]
            self.ignore_last_tstep = True
        else:
            raise ValueError(f'ntim differs: {ntims}')
        
        # list of all species
        varlists = [
                [att['VAR-LIST'][_*16:(_+1)*16].strip() for _ in range(len(att['VAR-LIST']) // 16)] 
                for att in atts]
        varlists2 = varlists[0].copy()
        for vl in varlists[1:]:
            #print(varlists2)
            varlists2 += [_ for _ in vl if _ not in varlists2]
        #print(varlists2)
        self.varlists2 = varlists2

        # check compatibility of species
        for v in varlists2:
            units = {self.fnames[ifile]:ff.variables[v].units for ifile, ff in enumerate(f) if v in ff.variables}
            shapes = {self.fnames[ifile]:ff.variables[v].shape for ifile, ff in enumerate(f) if v in ff.variables}

            units0 = next(iter(units.values())) 
            shape0 = next(iter(shapes.values())) 

            if not all(u == units0 for u in units.values()):
                unique_units = list(set(units.values()))
                units_and_files = {uu: [fn for  fn,un in units.items() if un==uu] for uu in unique_units}

                raise ValueError('Inconsitent Units: {units_and_files}')

            if not all(shp == shape0 for shp in shapes.values()):
                unique_shapes = list(set(shapes.values()))
                shapes_and_files = {shp: [fn for  fn,un in shapes.items() if shp==us] for us in unique_shapes}

                raise ValueError('Inconsitent Shapes: {shapes_and_files}')
            print(f'{v}: units = {units0}, shape = {shape0}')



        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

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
    p.add_argument('-s', '--summary', help='summary csv file name')

    args = p.parse_args()

    ptm = MrgUam(args.filenames, args.outname, args.summary)

if __name__ == '__main__':
    main()
