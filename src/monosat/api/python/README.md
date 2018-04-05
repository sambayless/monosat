# Monosat Python Interface
This package defines both a ctypes and (optional, but faster) cython API for MonoSAT
To install, first compile monosat, with
'''
cd ../../../../; cmake .; make
'''

Then install the python api with just ctypes support:
'''
$ sudo python3 setup.py install -f
'''

If cython is available on your system (eg, with pip3 install cython), as well an appropriate compiler, 
and 'use_cython=True' is set in setup.py, then this should compile and install the cython API.
Otherwise, it will install the (roughly 30% slower) ctypes api.
You can force the use of the ctypes or cython api, by manually setting the use_cython variable in either setup.py or monosat/monosat_c.py  
