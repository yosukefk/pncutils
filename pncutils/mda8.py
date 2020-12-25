import PseudoNetCDF as pnc

import numpy as np


class Mda8:

    def __init__(self, fnames, oname=None):
        """Moveing 8hr average

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

        self.buf = np.array(pnc.pncopen(fn).variables['O3'])

        i = 0
        for fn in self.fnames[1:]:
            # read the next file
            f = pnc.pncopen(fn)
            try:
                print(f.SDATE, i)
            except AttributeError:
                print(fn)
                raise

            o3 = np.array(f.variables['O3'])

            tmp = np.vstack([self.buf, o3])

            # get moving ave
            a8 = moving_average(tmp[:(24 + 7)])
            assert (a8.shape[0] == 24)
            a8 = a8.max(axis=0)

            # get rid of values at boudary cells
            a8[..., 0, :] = np.nan
            a8[..., -1, :] = np.nan
            a8[..., :, 0] = np.nan
            a8[..., :, -1] = np.nan

            # write
            self.fo.variables['MDA8O3'][i, ...] = a8[...]

            self.fo.updatetflag()

            self.buf = o3
            i += 1

    def _mkheader(self):

        f0 = pnc.pncopen(self.fnames[0])
        atts = f0.getncatts()

        # create new dataset
        fo = pnc.cmaqfiles.ioapi_base()

        # copy the dimensions
        for k, v in f0.dimensions.items():
            if k == 'TSTEP':
                n = len(self.fnames) - 1
            elif k == 'VAR':
                n = 1
            else:
                n = v.size
            fo.createDimension(v.name, n)
        fo.dimensions['TSTEP'].setunlimited(True)

        # set global attr
        atts['NVARS'] = 1
        atts['VAR-LIST'] = 'MDA8O3'.ljust(16)
        atts['TSTEP'] = 240000
        fo.setncatts(atts)

        # make variables
        fo.updatetflag()
        for i, nm in enumerate(['X', 'Y', 'longitude', 'latitude',
                                'O3']):
            v0 = f0.variables[nm]
            if i >= 4:
                prefix = 'MDA8'
            else:
                prefix = ''

            v = fo.createVariable(prefix + v0.name, v0.dtype.kind, v0.dimensions)
            a = v0.__dict__.copy()
            a['long_name'] = prefix + a['long_name']
            a['var_desc'] = prefix + ' ' + a['var_desc']
            v.setncatts(a)
            if i < 4:
                v[...] = v0[...]

        return fo


def _moving_average(x, w=8):
    return np.convolve(x, np.ones(w), 'valid') / w


def moving_average(x, w=8):
    """moving average on the first axis

    :param x: numpy float array
    :param w: (optional) window size, default 8
    :return: numpy float array, first dim shorter py w-1
    """
    n = x.shape[0]
    m = n - w + 1
    o = np.zeros((m,) + x.shape[1:])
    for i in range(n):
        for j in range(w):
            if i - j < 0 or i - j >= m: continue
            o[(i - j), ...] += x[i, ...]
    return o / w


def tester():
    o3 = Mda8(
        """
camx65ss_cb6r4hCF.20160101.rh.bc16_16s1.v1LN_v1a.2016_wrf381_p2KFsn_i2KFsn.avrg.grd01.nc
camx65ss_cb6r4hCF.20160102.rh.bc16_16s1.v1LN_v1a.2016_wrf381_p2KFsn_i2KFsn.avrg.grd01.nc
camx65ss_cb6r4hCF.20160103.rh.bc16_16s1.v1LN_v1a.2016_wrf381_p2KFsn_i2KFsn.avrg.grd01.nc
""".strip().split()
        , 'ooo.nc')
    return o3


# o3 = tester()

def main():
    import sys
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('-o', '--out', dest='oname', help='output file name')
    # opts, args = getopt.getopt(sys.argv[1:], 'o:')
    opts, args = p.parse_args()

    Mda8(args, opts.oname)


if __name__ == '__main__':
    main()
