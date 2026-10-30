"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
import os, re

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

VERSION_FILE = os.path.join(here, "p2rank", "version.py")

def get_version():
    with open(VERSION_FILE, "r") as f:
        content = f.read()
    match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', content, re.M)
    if not match:
        raise RuntimeError(f"Unable to find version string in {VERSION_FILE}")
    return match.group(1)

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-chem-p2rank',  # Required
    version=get_version(),  # Required
    description='Scipion plugin in order to use the p2rank software',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-chem/scipion-chem-p2rank',  # Optional
    author='Daniel Del Hoyo',  # Optional
    author_email='ddelhoyo@cnb.csic.es',  # Optional
    keywords='scipion p2rank scipion-3.0 cheminformatics',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={'pyworkflow.plugin': 'p2rank = p2rank'},
    package_data={  # Optional
       'p2rank': ['p2rank_logo.png', 'protocols.conf'],
    }
)
