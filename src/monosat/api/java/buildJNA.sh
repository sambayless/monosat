#!/bin/sh
#Use to build no-frills Java interface to MonoSAT's C API
#tested with jnaerator 0.9.4 and JNA 3.4
java -jar jnaerator.jar -DJNA=1  -direct -nocpp  -library  monosat Monosat.h
