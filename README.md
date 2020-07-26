# MonoSAT
[![Build Status](https://travis-ci.org/sambayless/monosat.svg?branch=master)](https://travis-ci.org/sambayless/monosat)


MonoSAT is a SAT Modulo Theory solver for *[monotonic theories]*, over Booleans and bitvectors. It supports a wide set of graph predicates (including reachability, shortest paths, maximum *s-t* flow, minimum spanning tree, and acyclicity constraints). 
MonoSAT supports reasoning about graphs that are either directed or undirected (including graphs that may have cycles). Edges may (optionally) have constant or Bitvector weights.
MonoSAT also has experimental support for constraints on finite state machines - see [FiniteStateMachines] for details.

MonoSAT can be used from the command line, or as a Python 3 library. See the installation instructions below; see also the [tutorial].

To see further examples of use cases for MonoSAT, and details on the (very simple) input file format that MonoSAT accepts, see  [FORMAT].

### Building
MonoSAT requires CMake (version 2.7 or higher).

From the root directory, build and install MonoSAT with:

```
$cmake .
$make
$sudo make install
```

MonoSAT requires C++11 support, zlib, and GMP >= 5.1.3. Tested on Ubuntu 14.04 and 16.04, with g++ 4.8.2 or higher, and with Clang 3.5. The Python library requires Python 3.3+.

If you get compilation errors along the lines of "error: could not convert ‘x’ from ‘__gmp_expr<__mpq_struct [1], __mpq_struct [1]>’ to ‘bool’", then you likely need to install a more recent version of GMP, with C++ support enabled.

If you build MonoSAT without using the provided cmake/makefiles, it is critically important to compile with `NDEBUG` set (*i.e.,* `-DNDEBUG`), as otherwise many very expensive debugging assertions will be enabled. 

#### Building on OSX with brew

You will need to first install GMP:

```
$brew install gmp
```

MonoSAT will, by default, compile statically linked binaries and dynamically linked libraries.
In most cases, you can then continue the cmake installation as above.

However, in some cases brew might not add the GMP libraries to the library path, in which case MonoSAT compilation may fail during linking.
You can check where brew installed GMP with:

```
$brew --prefix gmp
```

If brew installed GMP to /opt/local/lib/gmp, then you may need to modify the MonoSAT compilation instructions as follows:

```
$cmake .
$DYLD_LIBRARY_PATH=/opt/local/lib LIBRARY_PATH=/opt/local/lib make
$sudo make install
```

#### Building on FreeBSD *(Only tested on FreeBSD 12)*

You will need to first install GMP, which you can do via pkg or ports:

```
$pkg install gmp
```

If you intend to include the Java library (see below), you'll also need a JDK and you'll need to set your JAVA_HOME environment variable:

```
$pkg install openjdk8
$export JAVA_HOME=/usr/local/openjdk8
```

### Installing the Python library

To install the Python library (system-wide) on your system's default Python version:
```
$cmake -DPYTHON=ON .
$make
$sudo make install
```

MonoSAT also has support for Cython bindings, which are about 30% faster but require you to have a wokring installation of cython:
```
$cmake -DPYTHON=ON -DCYTHON=ON .
$make
$sudo make install
```

To install MonoSAT in an alternative python version (or in a virtual env or a non-standard location), first build MonoSAT as normal, then cd into 'src/monosat/api/python', and then use `setup.py` to install the Python library manually, eg:
```
$cd src/monosat/api/python
$sudo python3.6 setup.py install -f
```

Above, `-f` ensures that if you have previously installed a version of the monosat library, it will be overwritten cleanly with the new version.

See the [tutorial] and [tutorial.py] for instructions on using the Python library.

### Compiling the Java Library

To compile the Java library, you need to build MonoSAT with Java bindings enabled.
You will also need an installed JDK, version 1.8 or higher:

```
$cmake -DJAVA=ON .
$make
```

This should generate ```monosat.jar``` in MonoSAT's root directory.
To use MonoSAT from Java, you will need to include the jar in your classpath.
You will also need to ensure that Java can find MonoSAT's dynamic library, for example:
```
java -Djava.library.path=path/to/libmonosat.so -cp path/to/monosat.jar mypacakge.MyMainClass
```

On OSX, you would instead use ```path/to/libmonosat.dylib```

See [Tutorial.java] for instructions on using the Java library.

### Command-Line Usage

The recommended way to use MonoSAT is as a library, via the Python or Java bindings.
However, it is also possible to call MonoSAT from the command line.
MonoSAT is based on [MiniSat 2][Minisat], and supports many of the same calling conventions:

```
$monosat [-witness|-witness-file=filename] input_file.gnf
```

Where input_file.gnf is a file in [GNF format][FORMAT] (a very simple extension of DIMACS CNF format to support graph, bitvector, and finite state machine predicates). Use `-witness` to print the solution (if one exists) to stdout, or `-witness-file` to save it to file.

MonoSAT includes a very large set of configuration options - most of which you should stay away from unless you know what you are doing or want to explore the internals of MonoSAT (also, some of those configuration options might lead to buggy behaviour). Two options that often have a large impact on performance are `-decide-theories` and `-conflict-min-cut`:

```
$monosat -decide-theories -conflict-min-cut input_file.gnf
```

The `-decide-theories` option will cause the solver to make heuristic decisions that try to satisfy the various SMT predicates, which will often lead to improved performance, but can be pathologically bad in some common cases, and so is disabled by default. `-conflict-min-cut` will cause the solver to use a much slower, but more aggressive, clause learning strategy for reachability predicates; it may be useful for small, dificult instances.

MonoSAT implements a generalization of the circuit routing heuristics described in [Routing Under Constraints](#nadelruc16); you can activate them using the '-ruc' command line option. Thse can greatly improve performance on instances that use multiple reachability constraints. See [`examples/python/routing/router.py`][router] for an example usage.

### Source Overview
MonoSAT is written in C++. Core SAT solver functionality is in the `core/` and `simp/` directories; in particular, note `core/Config.cpp`, which is a central listing of all the configuration options available to MonoSAT. 

The graph and finite state machine theory solvers can be found in `graph/` and `fsm/`. Many of the graph algorithsms used by MonoSAT are collected in  `dgl/` (for 'Dynamic Graph Library').

`dgl/` incldudes C++ implementations of several dynamic graph algorithms (as well as some more common graph algorithms), and is well-optimized for medium sized (<20,000 nodes, < 100,000 edges), sparse graphs. The algorithms in dgl are designed for the case where the set of *possible* edges (and nodes) is fixed and known in advance (or only changes infrequently), and from that fixed set of possible edges many subsets of edges will subsequently be selected to be included in or excluded from the graph. 'dgl' supports templated edge weights and edge capacities, and has been tested successfully with integers, floats, and GMP arbitrary precision rationals.

##### Algirthms implemented in 'dgl/' include:
* Reachability/Shortest Path
    * [Ramalingam-Reps](#ramalingam1996incremental) dynamic single-source shortest paths algorithm (with improvements described in [Buriol et al. 2008](#buriol2008speeding))
    * [Thorup](#thorup2000near)'s dynamic connectivity algorithm (for undirected graphs)
    * [Dijkstra](#dijkstra1959note)'s algorithm
    * DFS/BFS
* Minimum Spanning Tree
    * [Spira and Pan](#spira1975finding)'s dynamic minimum spanning tree algorithm
    * [Kruskal](#kruskal)'s algorithm
    * [Prim](#prims)'s algorithm
* Maximum *s-t* Flow
    * [Kohli and Torr](#kohlitorr)'s dynamic variant of Kolmogorov and Boykov's maximum flow algorithm.
    * [Edmonds-Karp](#edmondskarp) algorithm (including a [dynamic variant](#dynamic_edmonds_karp))
    * [Dinitz](#dinitz)'s algorithm (including the [dynamic tree](#dynamic_tree) variant)
* And many supporting algorithms and data structures, including:
    * [Pearce and Kelly's `PK'](#pktopo) algorithm for dynamic topological sort in DAGs
    * [Euler Tree](http://en.wikipedia.org/wiki/Euler_tour_technique)
    * [Link Cut Tree](http://en.wikipedia.org/?title=Link/cut_tree)
    * [Tarjan's SCC](https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm)
    * Splay Tree
    * Disjoint Set


### Licensing
The majority of MonoSAT is released under the [MIT license] (as documented in individual source files). 
However, by default MonoSAT links some GPLv2 sources (found in ```src/monosat/dgl/alg/dyncut/```).
If built with these sources, the resulting binary is licensed as a whole under the [GPLv2].

MonoSAT can be built without including any GPL licensed sources, in which case it retains the MIT license.
To build MonoSAT without using GPL sources, use:
'''
$cmake -DGPL=OFF
''' 


### Acknowledgements

MonoSAT was made possible by the use of several open-source projects, including the afore-mentioned [MiniSat], as well as a high-performance [dynamic maximum-flow algorithm] by Pushmeet Kohli and Philip Torr, Emil Stefanov's implementation of [Disjoint Sets], and a [Link-Cut Tree] implementation by Daniel Sleator.

[monotonic theories]: http://www.cs.ubc.ca/labs/isd/Projects/monosat/smmt.pdf
[FORMAT]: FORMAT.md
[tutorial]: TUTORIAL.md
[FiniteStateMachines]: FiniteStateMachines.md
[tutorial.py]: examples/python/tutorial.py
[Tutorial.java]: examples/java/Tutorial.java
[router]: examples/python/routing/router.py
[MiniSat]:http://minisat.se/

[MIT license]: http://opensource.org/licenses/MIT
[GPLv2]: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
[dynamic maximum-flow algorithm]:http://research.microsoft.com/en-us/um/people/pkohli/code/rrr.txt
[Link-Cut Tree]: http://codeforces.com/contest/117/submission/860934
[Disjoint Sets]: http://web.rememberingemil.org/Projects/DisjointSets.aspx.html

### Publications using MonoSAT
* [Bayless, S., Bayless, N., Hoos, H. H., Hu, A.J. "SAT Modulo Monotonic Theories", Proceedings of the 29th AAAI Conference on Artificial Intelligence (2015)](http://www.cs.ubc.ca/labs/isd/Projects/monosat/smmt.pdf)
* [Klenze, T., and Bayless, S., and Hu, A.J. "Fast, Flexible, and Minimal CTL synthesis via SMT", International Conference on Computer Aided Verification (2016)](http://www.cs.ubc.ca/labs/isd/Projects/monosat/monosat_ctl_cav2016.pdf)
* <a name="nadelruc16">[Nadel, A. "Routing under Constraints." Formal Methods in Computer-Aided Design FMCAD (2016)](http://dl.acm.org/citation.cfm?id=3077653)</a>
* [Bayless, S., Hoos, H. H., and Hu, A.J. "Scalable, high-quality, SAT-based Multi-Layer Escape Routing", Proceedings of the 35th International Conference on Computer-Aided Design (2016)](http://www.cs.ubc.ca/labs/isd/Projects/monosat/monosat_bga_iccad2016.pdf)
* [Bayless, S. "SAT Modulo Monotonic Theories." PhD thesis, University of British Columbia (2017)](http://www.cs.ubc.ca/labs/isd/Projects/monosat/sam_bayless_thesis_2017.pdf)
* [Bayless, S., Kodirov, N., Beschastnikh, N., Hoos, H. H., Hu, A.J. "Scalable Constraint-Based Virtual Data Center Allocation", International Joint Conference on Artificial Intelligence (2017)](http://www.cs.ubc.ca/labs/isd/Projects/monosat/netsolver_ijcai2017.pdf)

If you would like your publication listed here, please contact [Sam Bayless](mailto:sbayless@cs.ubc.ca).

If you want to cite MonoSAT, please cite our 2015 AAAI paper:
```
@article{monosat2015,
  author	= {Sam Bayless and Noah Bayless and Holger H. Hoos and Alan J. Hu},
  title		= {{SAT Modulo Monotonic Theories}},
  booktitle	= {Proceedings of the 29th AAAI Conference on Artificial Intelligence},
  year		= {2015}
}
```

### Further References
* <a name="buriol2008speeding">[Buriol, Luciana S., Mauricio GC Resende, and Mikkel Thorup. "Speeding up dynamic shortest-path algorithms." INFORMS Journal on Computing 20.2 (2008): 191-204.](http://dx.doi.org/10.1287/ijoc.1070.0231)</a>
* <a name="dijkstra1959note">[Dijkstra, Edsger W. "A note on two problems in connexion with graphs." Numerische mathematik 1.1 (1959): 269-271](http://dx.doi.org/10.1007%2FBF01386390)</a>
* <a name="dinitz">[Dinitz, Y. "Algorithm for solution of a problem of maximum flow in a network with power estimation". Doklady Akademii nauk SSSR 11: 1277–1280  (1970)](http://www.cs.bgu.ac.il/~dinitz/D70.pdf)</a>
* <a name="edmondskarp">[Edmonds, Jack, and Richard M. Karp. "Theoretical improvements in algorithmic efficiency for network flow problems." Journal of the ACM (JACM) 19.2 (1972): 248-264](http://dx.doi.org/10.1145%2F321694.321699)</a>
* <a name="kohli2005efficiently">[Kohli, Pushmeet, and Philip HS Torr. "Efficiently solving dynamic markov random fields using graph cuts." Computer Vision, 2005. ICCV 2005. Tenth IEEE International Conference on. Vol. 2. IEEE, 2005](http://dx.doi.org/10.1109/ICCV.2005.81)</a>
* <a name="dynamic_edmonds_karp">[Korduban, D. "Incremental Maximum Flow in Dynamic graphs." Theoretical Computer Science Stack Exchange. http://cstheory.stackexchange.com/q/10186, 2012](http://cstheory.stackexchange.com/a/10186)</a>
* <a name="kruskal">[Kruskal, Joseph B. "On the shortest spanning subtree of a graph and the traveling salesman problem." Proceedings of the American Mathematical society 7.1 (1956): 48-50](http://dx.doi.org/10.1090%2FS0002-9939-1956-0078686-7)</a>
* <a name="prim">[Prim, Robert Clay. "Shortest connection networks and some generalizations." Bell system technical journal 36.6 (1957): 1389-140](http://dx.doi.org/10.1002/j.1538-7305.1957.tb01515.x)</a>
* <a name="ramalingam1996incremental"> [Ramalingam, Ganesan, and Thomas Reps. "An incremental algorithm for a generalization of the shortest-path problem." Journal of Algorithms 21.2 (1996): 267-305.](http://dx.doi.org/10.1006/jagm.1996.0046)</a>
* <a name="dynamic_tree">[Sleator, Daniel D., and Robert Endre Tarjan. "A data structure for dynamic trees." Proceedings of the thirteenth annual ACM symposium on Theory of computing. ACM, 1981.](http://dx.doi.org/10.1145/800076.802464)</a>
* <a name="spira1975finding">[Spira, Philip M., and A. Pan. "On finding and updating spanning trees and shortest paths." SIAM Journal on Computing 4.3 (1975): 375-380.](http://dx.doi.org/10.1137/0204032)</a>
* <a name="thorup2000near">[Thorup, Mikkel. "Near-optimal fully-dynamic graph connectivity." Proceedings of the thirty-second annual ACM symposium on Theory of computing. ACM, 2000.](http://dx.doi.org/10.1145/335305.335345)</a>
* <a name="pktopo">[Pearce, David J., and Kelly, Paul H. J. "A Dynamic Topological Sort Algorithm for Directed Acyclic Graphs", 2006](http://dx.doi.org/10.1145/1187436.1210590)</a>
