#!/usr/bin/env python3
from __future__ import print_function

import os.path
import platform
import shutil
import sys
try:
    from setuptools import setup
except ImportError as e:
    from distutils.core import setup


if sys.version_info[0] < 3:
    sys.exit('Sorry, Python < 3 is not supported')

# Set to False to disable compiling cython modules, set to True to enable cython

use_cython = False

# let cmake configure whether we use cython or not
# this string will be replaced by cmake with a string literal in that case
cmake_use_cython = '${USE_CYTHON}'

if cmake_use_cython.startswith("$"):
    pass  # cmake did not configure this file
else:
    use_cython = (cmake_use_cython == "True")

# allow the user to set whether cython is used using an environment variable
if "MONOSAT_CYTHON" in os.environ:
    use_cython = str(os.environ["MONOSAT_CYTHON"]) == "1"

# allow cmake to configure the package directory
package_dir = '${PACKAGE_DIR}'
if package_dir.startswith("$"):
    package_dir = '.'

library_dir = "${CMAKE_BINARY_DIR}"
if library_dir.startswith("$"):
    library_dir = "../../../../"

monosat_path = "${CMAKE_SOURCE_DIR}/src"
if monosat_path.startswith("$"):
    monosat_path = "../../../../src/"

if use_cython:
    print("Attempting Cython installation")
    # attempt to load the cython modules
    try:
        from distutils.extension import Extension
        from Cython.Build import cythonize
        from Cython.Distutils import build_ext
        from distutils.command.sdist import sdist as _sdist
    except:
        print("Could not load cython modules, falling back on ctypes")
        use_cython = False

if platform.system() == "Darwin":
    sharedlib = 'libmonosat.dylib'
elif platform.system() != "Windows":
    sharedlib = 'libmonosat.so'
else:
    sharedlib = 'libmonosat.dll'

orig_lib = library_dir + "/" + sharedlib
copy_lib = package_dir + "/monosat/" + sharedlib
if os.path.exists(orig_lib):
    # only copy the library if it hasn't already been copied (this facilitates separate build/install steps)
    if not os.path.exists(copy_lib) or os.path.getmtime(orig_lib) > os.path.getmtime(copy_lib):
        shutil.copy2(orig_lib, package_dir + "/monosat/")

if not os.path.exists(package_dir + "/monosat/" + sharedlib):
    print("Warning: could not find %s. See README for instructions on compiling the library, the re-install" % (
        sharedlib), file=sys.stderr)

if use_cython:

    # build the cython interface to monosat
    cmdclass = {}
    cmdclass.update({'build_ext': build_ext})
    setup(
        version='1.6',
        python_requires='>3.0.0',
        description='MonoSAT Cython Interface',
        author='Sam Bayless',
        author_email='sbayless@cs.ubc.ca',
        url='http://www.cs.ubc.ca/labs/isd/projects/monosat/',
        cmdclass=cmdclass,
        runtime_library_dirs=['./', package_dir + "/"],
        ext_modules=cythonize([Extension("monosat.monosat_p", [package_dir + "/monosat/monosat_p.pyx"],
                                         include_dirs=[".", package_dir, package_dir + "/monosat", monosat_path],
                                         libraries=["monosat"],
                                         language="c", extra_compile_args=["-DNDEBUG", "-O3"]
                                         )], include_path=[package_dir, package_dir + "/monosat"], gdb_debug=True),
        install_requires=['cython'],
        packages=['monosat'],
        package_data={'monosat': [sharedlib]},
        package_dir={'': package_dir},

    )
else:
    setup(name='monosat',
          version='1.6',
          python_requires='>3.0.0',
          description='MonoSAT Python Interface',
          author='Sam Bayless',
          author_email='sbayless@cs.ubc.ca',
          url='http://www.cs.ubc.ca/labs/isd/projects/monosat/',
          packages=['monosat'],
          package_data={'monosat': [sharedlib]},
          package_dir={'': package_dir},
          )
