import sys
from os.path import abspath, split, join
from setuptools import setup, find_packages

src_path = split(abspath(sys.argv[0]))[0]
src_path = join(src_path, 'mappgene')

# ------------------------------------------------------------------------------
setup(
    name='mappgene',
    description='Genomic sequence analysis for HPC',
    version='0.1.2',
    python_requires='>3.7.0',
    author='Joseph Moon',
    author_email='jmoon1506@gmail.com',
    entry_points={
        'console_scripts': [
            'mappgene = mappgene.cli:main'
        ]
    },
    packages=find_packages(),
    include_package_data=True,
    package_data={'mappgene' : ['data/*'] },
    install_requires=['numpy>=1.21.2', 'parsl>=1.1.0', 'pytest>=6.2.4'],
)

# ------------------------------------------------------------------------------
