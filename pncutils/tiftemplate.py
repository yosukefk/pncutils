'''takes a camx ncf as location template and create Gtiff template'''

import PseudoNetCDF as pnc
try:
    import gdal, osr
except ImportError:
    from osgeo import gdal, osr
import tempfile


class TifTemplate:
    '''takes a camx ncf as location template and create Gtiff template'''

    def __init__(self, fname):
        '''camx file as geo template to make your tif file

        :param fname: tempalte camx nc file (read only)


        EXAMPLE USAGE:
        probably METHOD1 makes more sense?


        # create template first
        tt = TifTemplate('my.nc')

        # either, (METHOD 1) specify pass and array to write, 
        tt.mktif('out.tif', my_2d_array)

        # or, (METHOD 2)get the raster definitions (template) and use it
        template = tt.mktif()
        o = tt.drv.CrateCopy('out.tif', template)
        o.GetRasterBand(1).WriteArray(my_wd_array)
        o.FlushCache()
        del o


            '''
        self.fname = fname

        self.drv = gdal.GetDriverByName('GTiff')

    def mktif(self, oname=None, arr=None, dtype=None):
        '''create tif file based on template

        :rtype: gdal.Dataset
        :param oname: output tif name
        :param arr: numpy 2d array to write to tif
        :param dtype: numpy dtype

        if oname is not provided, you can use returned gdal.Dataset to create tif by driver.CreateCopy()
        
        '''

        if oname is None:
            of = tempfile.NamedTemporaryFile()
            oname = of.name
        oname = str(oname)
        self.dtype = dtype

        ds = pnc.pncopen(self.fname)
        #nx, ny = [int(_) for _ in (ds.NCOLS, ds.NROWS)]
        nx, ny = [int(getattr(ds, _)) for _ in ('NCOLS', 'NROWS')]
        # pnc returns projection with false easting/norting.  sor gtrans shouldnt include xo yo info
        #gtrans = [float(_) for _ in (ds.XORIG, ds.XCELL, 0, ds.YORIG + ny * ds.YCELL, 0, - ds.YCELL, )]
        gtrans = [float(_) for _ in (0, ds.XCELL, 0, + ny * ds.YCELL, 0, - ds.YCELL, )]

        srs = osr.SpatialReference()
        srs.ImportFromProj4(ds.getproj().to_proj4())

        if not self.dtype is None:
            warnings.warn('dtype stuck with float32 for now')

        ds = self.drv.Create(oname, xsize=nx, ysize=ny, bands=1, eType=gdal.GDT_Float32 )

        ds.SetGeoTransform(gtrans)
        ds.SetProjection(srs.ExportToWkt())

        if not arr is None:
            ds.GetRasterBand(1).WriteArray(arr)
        ds.FlushCache()
        return ds

def tester():
    fname = 'camx_cb6_ei_lo_loptus.2019_day.txo3.bc19.negu_2019_v2.txs_4km.nc'
    return TifTemplate(fname).mktif()

if __name__ == '__main__':
    help(TifTemplate)

    
    




