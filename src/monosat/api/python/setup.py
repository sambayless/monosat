#!/usr/bin/env python3
from __future__ import print_function

import os.path
import platform
import shutil
import sys
from distutils.core import setup
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

setup(name='monosat',
      version='1.4',
      description='MonoSAT Python Interface',
      author='Sam Bayless',
      author_email='sbayless@cs.ubc.ca',
      url='http://www.cs.ubc.ca/labs/isd/projects/monosat/',
      packages=['monosat'],
      package_data={'monosat': [sharedlib]},
      )
