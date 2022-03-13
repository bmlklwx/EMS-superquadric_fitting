import setuptools
setuptools.setup(name='EMS',
version='0.1',
description='EMS: a package for probabilistic recovery of superquadrics from point clouds',
url='#',
author='Weixiao Liu',
install_requires=[
    'numpy',
    'scipy',
    'plyfile',
    'mayavi',
    'numba'
    ],
author_email='wliu72@jhu.edu',
packages=setuptools.find_packages(),
zip_safe=False)