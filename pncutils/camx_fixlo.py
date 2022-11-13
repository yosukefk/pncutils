# as of 2022-01, EPS3 generates gridded emission file which has dim LAY = 1 with VGLVLS=[0], which causes trouble
# this script checks that;s that case and modeify VGLVLS to [0,2]
import PseudoNetCDF as pnc
import subprocess

def chk(fname, oname):
    f = pnc.pncopen(fname)
    atts = f.getncatts()
    if atts['NLAYS'] == 1 and ('__len__' not in dir(atts['VGLVLS'])) :
        subprocess.run(shlex.split(f'ncatted -a VGLVLS,global,o,d,"0.0,0.0" {fname} {oname}'))

        



if __name__ == '__main__':
    import sys
    import subprocess
    import shlex
    fname = sys.argv[1]
    # TODO check 
    #subprocess.run(shlex.split("ncdump -h test.nc | sed -n '/:NTHIK = /s/^[^=]*= \([^;]\+\);/\1/p'"))
    #subprocss.run(shlex.split(f'ncatted -a VGLVLS,global,o,d,"0,2" {fname} {fname[:-3]}.fixed.nc')
    try:
        oname = sys.argv[2]
    except IndexError:
        oname = 'fixlo.nc'
    chk(fname, oname)

