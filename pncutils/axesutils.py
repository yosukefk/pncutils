'''miscelaneous unit conv'''

from pyproj import CRS, Transformer
import datetime

crs_emep = CRS.from_proj4('+proj=longlat +R=6370000 +no_defs')
crs_tceq = CRS.from_proj4('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs')
trns = Transformer.from_crs(crs_tceq, crs_emep, always_xy=True)

def ij_to_lnlt(i0,j0, ds):
    """python 2-d index to lat/lon

    param i0: 0-based easting index (netcdf loads data with this indexing)
    param j0: 0-based northing index
    param ds: ioapi file (to read the grid structure info)

    returns tuple i1, j1, x, y, ln, lt: 
        i1, j1: 1-based easting/northin index
        x, y: projected coordinate in meters
        ln, lt: gegraphical coodinate, assuming it uses TCEQ sphere/projection (same for EPA)
    """
    i1 = i0 + 1
    j1 = j0 + 1 
    
    x = ds.XORIG + ds.XCELL * (i0 + .5)
    y = ds.YORIG + ds.YCELL * (j0 + .5)
    lt, ln = trns.transform(x, y)
    return i1, j1, x, y, ln, lt


def tflag_to_time(tflag):
    """IOAPI timestamp to python timestamp

    :param tf: YYYYJJJ and HHMMSS
    :returns: python datetime
    """
    t = ( datetime.datetime(tflag[0] // 1000 , 1, 1) 
                + datetime.timedelta(days = tflag[0] % 1000 - 1) 
                + datetime.timedelta(seconds = (tflag[1] // 10000) * 3600 + ((tflag[1] % 10000) // 100) * 60 + (tflag[1] % 100))
                )
    return t


