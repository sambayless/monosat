# JNI-based Java API for Monosat
##Contents
* monosat/Solver.java: Object-Oriented SAT solver interface to Monosat
* monosat/Lit.java: Represents Boolean literal in formulas
* monosat/LBool.java: Represents potentially unassigned Boolean values (True, False, or Undef)
* monosat/Compare.java: Enum containing basic comparison types (gt, lt, eq, etc)
* monosat/Graph.java: Object-Oriented API for Monosat's graph theory
* monosat/BitVector.java: Object-Oriented API for Monosat's BitVector theory
* monosat/Logic.java: Static accessor methods for common logical operations/constructions
* monosat/MonosatJNI: Low-level, no-frills interface to Monosat's C API.

## Installation
By default, the shared library does not compile with JNI bindings.
To build with JNI support (requires a JDK installation, tested with OpenJDK 1.8):
```
$ cd <root folder of monosat>
$ cmake  -DJAVA=ON
$ make 
```

This will compile JNI bindings into the shared library (libmonosat.so/dylib).
It will also produce monosat.jar in the root folder, containing the java API.

## Usage
In order to use the library, you will need to ensure that monosat.jar
is on the classpath. You will also need to ensure that libmonosat.so is in Java's
native library path, eg, java -Djava.library.path=/path/to/libmonosat.so

See examples/java/Tutorial.java for an introduction to using the library.

