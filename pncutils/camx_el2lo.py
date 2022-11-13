# convert elevated emission to lo local format

import PseudoNetCDF as pnc

class El2lo:

    def __init__(self, fname, loname, oname=None):
        """camx elevated to lolevel

        :param fname: input elevated file name
        :param loname: reference gridded file name
        :param oname: (optional) output file name
        """

        self.fname = fname
        self.loname = loname

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
        # attributes from existing files
        f0 = pnc.pncopen(self.fname)
        atts0 = f0.getncatts()
        f1 = pnc.pncopen(self.loname)
        atts1 = f1.getncatts()

        atts = self.fo.getncatts()

        x = f0.variables['xcoord']
        y = f0.variables['ycoord']
        idx = ((x[:]  - atts['XORIG'] ) / atts['XCELL']).astype(int)
        jdx = ((y[:]  - atts['YORIG'] ) / atts['YCELL']).astype(int)

        oob = (idx<0) | (idx>=atts['NCOLS']) | (jdx<0) | (jdx >= atts['NROWS'])
        print(sum(oob), len(oob))
        self.idx = idx
        self.jdx = jdx
        self.oob = oob

        species = [atts['VAR-LIST'][_*16:(_+1)*16].strip() for _ in range(len(atts['VAR-LIST']) // 16)]
        for s in species:
            print(s)
            try:
                v0 = f0.variables[s]
            except KeyError:
                continue
            for c in range(v0.shape[-1]):
                if oob[c]: continue
                self.fo.variables[s][:, 0, jdx[c], idx[c]] += v0[:, c]




        

    def _mkheader(self):

        # attributes from existing files
        f0 = pnc.pncopen(self.fname)
        atts0 = f0.getncatts()
        f1 = pnc.pncopen(self.loname)
        atts1 = f1.getncatts()

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
        print(1, atts)
        atts.update({k:v for k,v in atts1.items() if k in 
            ('NCOLS', 'NROWS', 'XORIG', 'YORIG', 'XCELL', 'YCELL', 
                    'CAMx_NAME', 'FILEDESC')})
        print(2, atts)
        fo.setncatts(atts)
        print(3, fo.getncatts())


        # make variables
        fo.updatetflag()
        #names = ['X', 'Y', 'TFLAG', 'ETFLAG', 'longitude', 'latitude'] 
        #names = ['X', 'Y',                    'longitude', 'latitude'] 
        names0 = ['X', 'Y', 'TFLAG', 'ETFLAG', 'longitude', 'latitude'] 

        names = names0 + [atts['VAR-LIST'][_*16:(_+1)*16].strip() for _ in range(len(atts['VAR-LIST']) // 16)]
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

            v = fo.createVariable(vv.name, vv.dtype.kind, dim)
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
        print(4, attso)
        return fo




def tester():

    fin = 'camx_cb6p_ei_el_neguus.2019_day.txo3.bc19.2019_pig_v2.nc'
    flo = 'camx_cb6_ei_lo_loptus.2019_day.txo3.bc19.negu_2019_v2.txs_4km.nc'
    fout = 'el2lo.nc'

    lo = El2lo(fin, flo, fout)
    return lo

lo = tester()



