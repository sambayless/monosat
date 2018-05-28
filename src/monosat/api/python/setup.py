#!/usr/bin/env python3
from __future__ import print_function

import os.path
import platform
import shutil
import sys
from distutils.core import setup

if sys.version_info[0] < 3:
    sys.exit('Sorry, Python < 3 is not supported')

monosat_path ="../../../../src/"
#Set to False to disable compiling cython modules, set to True to enable cython
use_cython=False

if use_cython:
    # attempt to load the cython modules
    try:
        from distutils.extension import Extension
        from Cython.Build import cythonize
        from Cython.Distutils import build_ext
        from distutils.command.sdist import sdist as _sdist
        use_cython=True
    except:
        print("Could not load cython modules, falling back on ctypes")
        use_cython=False

if platform.system() == "Darwin":
    sharedlib='libmonosat.dylib'
    copy_lib = "monosat/libmonosat.dylib"
    orig_lib = "../../../../libmonosat.dylib"

    if not os.path.exists(orig_lib) and os.path.exists("../../../../SharedLibrary/libmonosat.dylib"):
        orig_lib = "../../../../SharedLibrary/libmonosat.dylib"
    if not os.path.exists(orig_lib) and os.path.exists("../../../../OSX_SharedLibrary/libmonosat.dylib"):
        orig_lib = "../../../../OSX_SharedLibrary/libmonosat.dylib"

    if os.path.exists(orig_lib):
        if not os.path.exists(copy_lib):
            print("Copying %s to monosat/libmonosat.dylib"%(copy_lib))
            shutil.copy2(orig_lib, "monosat")
        if  os.path.getmtime(orig_lib) > os.path.getmtime(copy_lib):
            print("Copying %s to monosat/libmonosat.dylib"%(orig_lib))
            shutil.copy2(orig_lib, "monosat")

elif platform.system() != "Windows":
    sharedlib='libmonosat.so'
    copy_lib = "monosat/libmonosat.so"
    orig_lib = "../../../../libmonosat.so"

    if not os.path.exists(orig_lib) and os.path.exists("../../../../SharedLibrary/libmonosat.so"):
        orig_lib = "../../../../SharedLibrary/libmonosat.so"
    if not os.path.exists(orig_lib) and os.path.exists("../../../../OSX_SharedLibrary/libmonosat.so"):
        orig_lib = "../../../../OSX_SharedLibrary/libmonosat.so"

    if os.path.exists(orig_lib):
        if not os.path.exists(copy_lib):
            print("Copying %s to monosat/libmonosat.so"%(copy_lib))
            shutil.copy2(orig_lib, "monosat")
        if  os.path.getmtime(orig_lib) > os.path.getmtime(copy_lib):
            print("Copying %s to monosat/libmonosat.so"%(orig_lib))
            shutil.copy2(orig_lib, "monosat")

else:
    sharedlib='libmonosat.dll'
    copy_lib = "monosat/libmonosat.dll"
    orig_lib = "../../../../libmonosat.dll"
    if not os.path.exists(orig_lib) and os.path.exists("../../../../Win64SharedLibrary/libmonosat.dll"):
        orig_lib = "../../../../Win64SharedLibrary/libmonosat.dll"
    if os.path.exists(orig_lib):
        if not os.path.exists(copy_lib):
            print("Copying %s to monosat/libmonosat.so"%(copy_lib))
            shutil.copy2(orig_lib, "monosat")
        if  os.path.getmtime(orig_lib) > os.path.getmtime(copy_lib):
            print("Copying %s to monosat/libmonosat.so"%(orig_lib))
            shutil.copy2(orig_lib, "monosat")

if not os.path.exists(copy_lib):
    print("Warning: could not find libmonosat.so or libmonosat.dll. See README for instructions on compiling the library, the re-install",file=sys.stderr)

if use_cython:

    #build the cython interface to monosat
    cmdclass = { }
    cmdclass.update({ 'build_ext': build_ext })
    setup(
        version='1.4',
        python_requires='>3.0.0',
        description='MonoSAT Cython Interface',
        author='Sam Bayless',
        author_email='sbayless@cs.ubc.ca',
        url='http://www.cs.ubc.ca/labs/isd/projects/monosat/',
        packages=['monosat'],
        package_data={'monosat': [sharedlib]},
        cmdclass = cmdclass,
        runtime_library_dirs=['./'],
        ext_modules = cythonize([Extension("monosat.monosat_p", ["monosat/monosat_p.pyx"],
                                           include_dirs=[".","monosat",monosat_path], libraries=["monosat"],
                                           language="c",extra_compile_args=["-DNDEBUG","-O3"]
                                           )],include_path = ['.',"monosat"], gdb_debug=True),
        install_requires=['cython']

    )
else:
    setup(name='monosat',
          version='1.4',
          python_requires='>3.0.0',
          description='MonoSAT Python Interface',
          author='Sam Bayless',
          author_email='sbayless@cs.ubc.ca',
          url='http://www.cs.ubc.ca/labs/isd/projects/monosat/',
          packages=['monosat'],
          package_data={'monosat': [sharedlib]},
          )
