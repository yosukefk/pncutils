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
        )
