from setuptools import setup, find_packages, Extension
import re

modcontents = open('lmtslr/__init__.py').read()
version = re.search(r"__version__ = '([^']*)'",modcontents).group(1)

setup(
    name = "SpectralLineReduction",
    version = version,
    description = "SpectralLineReduction: Package for spectral line reduction for LMT",
    author = "Gopal Narayanan",
    author_email = "gopal@astro.umass.edu",
    packages = find_packages(),
    scripts = ['grid_data.py', 'multiprocess_otf_map.py', 'process_map.py',
               'process_otf_map.py', 'process_ps.py', 'view_cube.py',
               'view_spec_file.py']
    )
