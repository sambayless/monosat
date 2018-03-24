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

If cython is available on your system (eg, with pip3 install cython), as well an appropriate compiler, this should compile and install the cython API.
If cython is not available, then it will install with the (roughly 30% slower) ctypes api.
You can also force the use of the ctypes api, by manually setting the use_cython variable in either setup.py or monosat/monosat_c.py to False  
