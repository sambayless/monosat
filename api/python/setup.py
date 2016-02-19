#!/usr/bin/env python3

from distutils.core import setup

import os.path
import shutil
import platform

if platform.system != "Windows":
    copy_lib = "monosat/libmonosat.so"
    orig_lib = "../../SharedLibrary/libmonosat.so"
    if not os.path.exists(copy_lib):
        shutil.copy2(orig_lib, "monosat")
    if os.path.getmtime(orig_lib) > os.path.getmtime(copy_lib):
        shutil.copy2(orig_lib, "monosat")

setup(name='monosat',
      version='1.0',
      description='MonoSAT Python Interface',
      author='Sam Bayless',
      author_email='sbayless@cs.ubc.ca',
      url='http://www.cs.ubc.ca/labs/isd/Projects/monosat/',
      packages=['monosat'],
      package_data={'monosat': ['libmonosat.so']},
     )
