# try using distribute or setuptools or distutils.
try:
    import distribute_setup
    distribute_setup.use_setuptools()
except ImportError:
    pass

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
import sys
import re

# parse version from package/module without importing or evaluating the code
with open('pyani/__init__.py') as fh:
    for line in fh:
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group('version')
            break

if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: pyani requires Python 3 " +
                     "or above...exiting.\n")
    sys.exit(1)

setup(
    name="pyani",
    version=version,
    author="Leighton Pritchard",
    author_email="leighton.pritchard@hutton.ac.uk",
    description=''.join(["pyani provides a package and script for " +
                         "calculation of genome-scale average nucleotide " +
                         "identity."]),
    license="MIT",
    keywords="genome bioinformatics sequence",
    platforms="Posix; MacOS X",
    url="http://widdowquinn.github.io/pyani/",  # project home page
    download_url="https://github.com/widdowquinn/pyani/releases",
    scripts=[os.path.join('bin', 'average_nucleotide_identity.py'),
             os.path.join('bin', 'genbank_get_genomes_by_taxon.py')],
    packages=['pyani'],
    package_data={'pyani': ['tests/test_JSpecies/*.tab']},
    include_package_date=True,
    install_requires=['biopython',
                      'matplotlib',
                      'pandas',
                      'scipy',
                      'seaborn'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )
