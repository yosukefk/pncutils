from setuptools import setup

setup(
        name='pncutils',
        version='1.0.0',
        description='tools using PseudoNetCDF',
        url='https://github.com/yosukefk/pncutils',
        license='MIT',
        packages=['pncutils'], 
        install_requires=[
            'PseudoNetCDF', 
            ],
        entry_points={
            'console_scripts': [
                'camx_el2lo=pncutils.camx_el2lo:main',
                'camx_elwindow=pncutils.camx_elwindow:main',
                'camx_ptsmrg=pncutils.camx_ptsmrg:main',
                'camx_mrguam=pncutils.camx_mrguam:main',
                'camx_window=pncutils.camx_window:main',
                'camx_fixlo=pncutils.camx_fixlo:main',
                ],
            },
        )
