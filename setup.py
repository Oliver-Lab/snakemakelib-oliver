# Modelled on snakemakelib setup script
# --------------------------------------------------
# Imports
# --------------------------------------------------

from __future__ import print_function

# stdlib
import os
from setuptools import setup
from os.path import realpath, dirname, relpath, join

# --------------------------------------------------
# globals and constants
# --------------------------------------------------

ROOT = dirname(realpath(__file__))

# --------------------------------------------------
# classes and functions
# --------------------------------------------------

package_data = []


def package_path(path, filters=()):
    if not os.path.exists(path):
        raise RuntimeError("packaging non-existent path: %s" % path)
    elif os.path.isfile(path):
        package_data.append(relpath(path, 'snakemakelib-oliver'))
    else:
        for path, dirs, files in os.walk(path):
            path = relpath(path, 'snakemakelib-oliver')
            for f in files:
                if not filters or f.endswith(filters):
                    package_data.append(join(path, f))

rule_suffixes = ('.rules', '.rule')
workflow_suffixes = ('.workflow')
                    
package_path(join(ROOT, 'bin'))
package_path(join(ROOT, 'snakemakelib_oliver', '_templates'))
package_path(join(ROOT, 'snakemakelib_oliver', 'tools'))
package_path(join(ROOT, 'snakemakelib_oliver', 'rules'), rule_suffixes)
scripts = ['bin/slurmSnake']

REQUIRES = [
]

try:
    # Hack for readthedocs
    if not 'readthedocs' in os.path.dirname(os.path.realpath(__file__)):
        REQUIRES.append('matplotlib>=1.4.0')
    else:
        print("readthedocs in path name; assuming we're building docs @readthedocs")
        REQUIRES.append('sphinx-bootstrap-theme')
except:
    pass
    

# https://pythonhosted.org/setuptools/setuptools.html
SETUP_REQUIRES = [
]    

# Adding github to setup:
DEPENDENCY_LINKS = [
]

# Integrating pytest with setuptools: see
# https://pytest.org/latest/goodpractises.html#integrating-with-distutils-python-setup-py-test
from distutils.core import setup, Command
# you can also import from setuptools

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

setup(
    name="snakemakelib-oliver",
    author="Justin Fear",
    author_email="justin.m.fear@gmail.com",
    description="Snakemake rule library built on top of the snakemakelib packages",
    license="MIT",
    url="http://github.com/Oliver-Lab/snakemakelib-oliver",
    scripts=scripts,
    packages=[],
    package_data={'snakemakelib_oliver': package_data},
    data_files=[('config', ['data/slurm_cluster_config.yaml'])],
    setup_requires=SETUP_REQUIRES,
    install_requires=REQUIRES,
    dependency_links=DEPENDENCY_LINKS,
)
