import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='EMS',
    version='0.0.1',
    description='EMS: a package for probabilistic recovery of superquadrics from point clouds',
    url='https://github.com/bmlklwx/EMS-probabilistic_superquadric_fitting.git',
    author='Weixiao Liu, Yuwei Wu, Sipu Ruan, Gregory Chirikjian',
    author_email='wliu72@jhu.edu',
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    install_requires=[
        'numpy',
        'scipy',
        'plyfile',
        'mayavi',
        'numba'
    ],
    classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
    ],
    
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires='>=3.6'
)
