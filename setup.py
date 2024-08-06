from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

extensions = [
    Extension(
        'dp4.ConfPrune',
        sources=['dp4/ConfPrune.pyx'],
        include_dirs=[numpy.get_include()]
    )
]

setup(
    name='dp4',
    version='1.0',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'dp4': ['data/*']
    },
    install_requires=[
        'numpy', 'scipy', 'networkx', 'rdkit', 'setuptools'
    ],
    entry_points={
        'console_scripts': [
            'pydp4 = dp4.main:main', 
        ],
    },
    author='Kristaps Ermanis, Alexander Howarth, Jonathan M. Goodman, Benji Rowlands',
    description='PyDP4 integrated workflow for MM, DFT GIAO calculations, and DP4 analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/jbr46/DP4',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    ext_modules=cythonize(extensions),
)
