# JNI-based Java API for Monosat
##Contents
* monosat/MonosatJNI: Low-level, no-frills interface to Monosat's C API. 

## Installation
By default, the shared library does not compile with JNI bindings.
To build with JNI support:
```
$ cd <root folder of monosat>
$ cmake  -DJAVA_SUPPORT=ON
$ make 
```



