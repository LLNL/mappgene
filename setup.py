import sys
from pathlib import Path
from setuptools import setup, find_packages

this_dir = Path(__file__).parent
src_path = this_dir / "mappgene"
long_desc = (this_dir / "README.md").read_text()

# ------------------------------------------------------------------------------
setup(
    name='mappgene',
    description='Genomic sequence analysis for high-performance computing',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    version='1.3.0',
    python_requires='>=3.7',
    author='Joseph Moon',
    maintainer='Aram Avila-Herrera',  # TODO: replace with mailing list?
    maintainer_email='avilaherrera1@llnl.gov',
    entry_points={
        'console_scripts': [
            'mappgene = mappgene.cli:main'
        ]
    },
    packages=find_packages(),
    include_package_data=True,
    package_data={'mappgene': ['data/*']},
    install_requires=['pip>=21.2.4', 'parsl>=1.1.0', 'pytest>=6.2.4'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: System :: Distributed Computing",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)

# ------------------------------------------------------------------------------
