import sys
from os.path import abspath, split, join
from setuptools import setup, find_packages

src_path = split(abspath(sys.argv[0]))[0]
src_path = join(src_path, 'mappgene')

# ------------------------------------------------------------------------------
setup(
    name='mappgene',
    description='Genomic sequence analysis for HPC',
    version='0.1.5',
    python_requires='>=3.7',
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
    install_requires=['numpy>=1.11.3', 'parsl>=0.9.0', 'pytest>=4.0'],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)

# ------------------------------------------------------------------------------
