cmake_minimum_required(VERSION 2.7)
project(monosat)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare  -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS")

set( CMAKE_VERBOSE_MAKEFILE on )

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   message(STATUS "Setting build type to 'Release' as none was specified.\nUse cmake -DCMAKE_BUILD_TYPE=Debug" to seelct Debug release.)
   set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
   # Set the possible values of build type for cmake-gui
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"  "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "build type is " ${CMAKE_BUILD_TYPE})
include_directories(./)

set(SOURCE_FILES
    amo/AMOParser.h
    amo/AMOTheory.h
    api/Circuit.h
    api/Monosat.cpp
    api/Monosat.h
    bv/BVParser.h
    bv/BVTheory.h
    bv/BVTheorySolver.h
    core/AssumptionParser.h
    core/Config.cpp
    core/Config.h
    core/Dimacs.h
    core/Optimize.cpp
    core/Optimize.h
    core/Remap.h
    core/Solver.cc
    core/Solver.h
    core/SolverTypes.h
    core/Theory.h
    core/TheorySolver.h
    dgl/alg/dyncut/block.h
    dgl/alg/dyncut/graph.h
    dgl/alg/AugmentedSplayTree.h
    dgl/alg/DisjointSets.cpp
    dgl/alg/DisjointSets.h
    dgl/alg/EulerTree.h
    dgl/alg/Heap.h
    dgl/alg/LinkCut.h
    dgl/alg/LinkCutCost.cpp
    dgl/alg/LinkCutCost.h
    dgl/alg/SearchTree.h
    dgl/alg/SplayTree.h
    dgl/alg/Treap.h
    dgl/alg/TreapCustom.h
    dgl/AcyclicFlow.h
    dgl/AllPairs.cpp
    dgl/AllPairs.h
    dgl/BFS.h
    dgl/Chokepoint.h
    dgl/ConnectedComponents.cpp
    dgl/ConnectedComponents.h
    dgl/Cycle.h
    dgl/DFS.h
    dgl/DFSCycle.h
    dgl/Dijkstra.h
    dgl/DijkstraAllPairs.h
    dgl/Dinics.h
    dgl/DinicsLinkCut.h
    dgl/DisjointSetConnectedComponents.h
    dgl/Distance.cpp
    dgl/Distance.h
    dgl/DynamicConnectivity.h
    dgl/DynamicConnectivityImpl.h
    dgl/DynamicGraph.h
    dgl/DynamicNodes.h
    dgl/EdmondsKarp.h
    dgl/EdmondsKarpAdj.h
    dgl/EdmondsKarpDynamic.h
    dgl/FloydWarshall.h
    dgl/KohliTorr.h
    dgl/Kruskal.h
    dgl/MaxFlow.h
    dgl/MinimumSpanningTree.h
    dgl/NaiveDynamicConnectivity.h
    dgl/PKTopologicalSort.h
    dgl/Prim.h
    dgl/RamalReps.h
    dgl/Reach.cpp
    dgl/Reach.h
    dgl/SpiraPan.h
    dgl/SteinerApprox.h
    dgl/SteinerTree.h
    dgl/TarjansSCC.h
    dgl/ThorupDynamicConnectivity.h
    fsm/alg/NFAAccept.h
    fsm/alg/NFAGenerate.h
    fsm/alg/NFALinearGeneratorAcceptor.h
    fsm/alg/NFATransduce.h
    fsm/alg/NFATypes.h
    fsm/alg/NP0LAcceptor.h
    fsm/alg/NFAAcceptor.h
    fsm/alg/Dawg.h
    fsm/DynamicFSM.h
    fsm/FSMAcceptDetector.cpp
    fsm/FSMAcceptDetector.h
    fsm/FSMDawgAcceptDetector.cpp
    fsm/FSMDawgAcceptDetector.h
    fsm/FSMDetector.h
    fsm/FSMGeneratesDetector.cpp
    fsm/FSMGeneratesDetector.h
    fsm/FSMGeneratorAcceptorDetector.cpp
    fsm/FSMGeneratorAcceptorDetector.h
    fsm/FSMParser.h
    fsm/FSMTheory.h
    fsm/FSMTransducesDetector.cpp
    fsm/FSMTransducesDetector.h
    fsm/LSystem.h
    fsm/LSystemParser.h
    fsm/LSystemTheory.h
    fsm/P0LAcceptDetector.cpp
    fsm/P0LAcceptDetector.h
    geometry/alg/Heap.h
    geometry/bounds/BoundingBox.h
    geometry/bounds/BoundingSphere.h
    geometry/bounds/BoundingVolume.h
    geometry/cevans/halfspaceContainer.h
    geometry/cevans/halfspaceD2.h
    geometry/cevans/halfspaceD3.h
    geometry/cevans/indextable.h
    geometry/cevans/mathlib.cpp
    geometry/cevans/mathlib.h
    geometry/cevans/partitionspace.h
    geometry/cevans/point.h
    geometry/cevans/quickhull2D.h
    geometry/cevans/quickhull3D.h
    geometry/cevans/typedefs.h
    geometry/cevans/typeop.h
    geometry/cevans/zero.h
    geometry/polypartition/polypartition.cpp
    geometry/polypartition/polypartition.h
    geometry/ConvexHull.h
    geometry/ConvexHullCollisionDetection.cpp
    geometry/ConvexHullCollisionDetector.h
    geometry/ConvexHullDetector.cpp
    geometry/ConvexHullDetector.h
    geometry/ConvexPolygon.cpp
    geometry/ConvexPolygon.h
    geometry/Delaunay.h
    geometry/DelaunayPolypartition.cpp
    geometry/DelaunayPolypartition.h
    geometry/GeometryDetector.h
    geometry/GeometryParser.h
    geometry/GeometryTheory.h
    geometry/GeometryTypes.h
    geometry/GridHeightmap.cpp
    geometry/GridHeightmap.h
    geometry/Heightmap.h
    geometry/HeightmapDetector.h
    geometry/Line.h
    geometry/LineSegment.h
    geometry/MonotoneConvexHull.h
    #/MonotoneDelaunay.cpp
    #geometry/MonotoneDelaunay.h
    geometry/PointSet.cpp
    geometry/PointSet.h
    geometry/Polygon.h
    geometry/PolygonSet.h
    geometry/QuickConvexHull.cpp
    geometry/QuickConvexHull.h
    geometry/Shape.h
    graph/AllPairsDetector.cpp
    graph/AllPairsDetector.h
    graph/ConnectedComponentsDetector.cpp
    graph/ConnectedComponentsDetector.h
    graph/CycleDetector.cpp
    graph/CycleDetector.h
    graph/Detector.h
    graph/DistanceDetector.cpp
    graph/DistanceDetector.h
    graph/GraphParser.h
    graph/GraphTheory.h
    graph/GraphTheoryTypes.h
    graph/MaxflowDetector.cpp
    graph/MaxflowDetector.h
    graph/MSTDetector.cpp
    graph/MSTDetector.h
    graph/ReachDetector.cpp
    graph/ReachDetector.h
    graph/SteinerDetector.cpp
    graph/SteinerDetector.h
    graph/WeightedDijkstra.h
    graph/WeightedDistanceDetector.cpp
    graph/WeightedDistanceDetector.h
    mtl/Alg.h
    mtl/Alloc.h
    mtl/Bitset.h
    mtl/Deque.h
    mtl/Heap.h
    mtl/IntMap.h
    mtl/IntTypes.h
    mtl/Map.h
    mtl/Queue.h
    mtl/Rnd.h
    mtl/Sort.h
    mtl/Vec.h
    mtl/XAlloc.h
    pb/PbParser.h
    pb/PbTheory.h
    simp/SimpSolver.cc
    simp/SimpSolver.h
    utils/Options.cc
    utils/Options.h
    utils/ParseUtils.h
    utils/System.cc
    utils/System.h
    Main.cc)

add_executable(monosat ${SOURCE_FILES})
set_target_properties(monosat PROPERTIES EXCLUDE_FROM_ALL 1)
target_link_libraries(monosat z)
target_link_libraries(monosat m)
target_link_libraries(monosat gmpxx)
target_link_libraries(monosat gmp)

if (${UNIX})
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    target_link_libraries(monosat rt)
endif (${UNIX})

add_executable(monosat_static ${SOURCE_FILES})
set_target_properties(monosat_static PROPERTIES OUTPUT_NAME monosat)
target_link_libraries (monosat_static "-static-libgcc -static-libstdc++")
target_link_libraries(monosat_static z.a)
target_link_libraries(monosat_static m.a)
target_link_libraries(monosat_static gmpxx.a)
target_link_libraries(monosat_static gmp.a)

if (${UNIX})
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    target_link_libraries(monosat_static rt) #not static
endif (${UNIX})

add_library(libmonosat SHARED ${SOURCE_FILES} )
set_target_properties(libmonosat PROPERTIES OUTPUT_NAME monosat)
target_link_libraries(libmonosat z)
target_link_libraries(libmonosat m)
target_link_libraries(libmonosat gmpxx)
target_link_libraries(libmonosat gmp)
if (${UNIX})
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    target_link_libraries(libmonosat rt)
endif (${UNIX})


add_library(libmonosat_static ${SOURCE_FILES})
set_target_properties(libmonosat_static PROPERTIES OUTPUT_NAME monosat)
target_link_libraries(libmonosat_static z.a)
target_link_libraries(libmonosat_static m.a)
target_link_libraries(libmonosat_static gmpxx.a)
target_link_libraries(libmonosat_static gmp.a)
if (${UNIX})
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    target_link_libraries(libmonosat_static rt)
endif (${UNIX})

install(TARGETS monosat_static libmonosat libmonosat_static
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        #INCLUDES DESTINATION include/monosat
        )
install(DIRECTORY . DESTINATION include/monosat
        FILES_MATCHING PATTERN "*.h"
        PATTERN api/python EXCLUDE
        PATTERN api/java EXCLUDE
        PATTERN tests EXCLUDE
        REGEX "^\\..*" EXCLUDE #exclude hidden files
        PATTERN CMakeFiles EXCLUDE
        PATTERN examples EXCLUDE
        #exclude old build files, in case they are still in the source directory
        PATTERN Debug EXCLUDE
        PATTERN DbgRelease EXCLUDE
        PATTERN Release EXCLUDE
        PATTERN Static EXCLUDE
        PATTERN ClangDebug EXCLUDE
        PATTERN ClangRelease EXCLUDE
        PATTERN OSX_ClangRelease EXCLUDE
        PATTERN OSX_SharedLibrary EXCLUDE
        PATTERN SharedLibrary EXCLUDE
        PATTERN SharedLibraryDbg EXCLUDE
        PATTERN OSX_SharedLibrary EXCLUDE
        PATTERN Win32 EXCLUDE
        PATTERN Win64 EXCLUDE
        PATTERN Win64SharedLibrary EXCLUDE)
