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
To build with JNI support:
```
$ cd <root folder of monosat>
$ cmake  -DJAVA=ON
$ make 
```



