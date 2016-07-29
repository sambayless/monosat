cmake_minimum_required(VERSION 2.7)
project(monosat)

set (PB_SUPPORT 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3  -DNO_GMP -std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare  -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS")

set( CMAKE_VERBOSE_MAKEFILE on )

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   message(STATUS "Setting build type to 'Release' as none was specified.\nUse cmake -DCMAKE_BUILD_TYPE=Debug to seelct Debug release.")
   set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
   # Set the possible values of build type for cmake-gui
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"  "MinSizeRel" "RelWithDebInfo" "Sanitize")
endif()

set(CMAKE_CONFIGURATION_TYPES "Release;Debug;MinSizeRel;RelWithDebInfo;Sanitize" CACHE STRING "" FORCE)

message("* Adding build types...")
SET(CMAKE_CXX_FLAGS_SANITIZE
        "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used by the C++ compiler during address sanetize  builds."
        FORCE )
SET(CMAKE_C_FLAGS_SANITIZE
        "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used by the C compiler during address sanetize builds."
        FORCE )
SET(CMAKE_EXE_LINKER_FLAGS_SANITIZE
        "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used for linking binaries during address sanetize builds."
        FORCE )
SET(CMAKE_SHARED_LINKER_FLAGS_SANITIZE
        "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used by the shared libraries linker during address sanetize builds."
        FORCE )
MARK_AS_ADVANCED(
        CMAKE_CXX_FLAGS_SANETIZE
        CMAKE_C_FLAGS_SANETIZE
        CMAKE_EXE_LINKER_FLAGS_SANETIZE
        CMAKE_SHARED_LINKER_FLAGS_SANETIZE )

include_directories(src/)

set(SOURCE_FILES 
    src/monosat/amo/AMOParser.h
    src/monosat/amo/AMOTheory.h
    src/monosat/api/Circuit.h
    src/monosat/api/Monosat.cpp
    src/monosat/api/Monosat.h
    src/monosat/bv/BVParser.h
    src/monosat/bv/BVTheory.h
    src/monosat/bv/BVTheorySolver.h
    src/monosat/core/AssumptionParser.h
    src/monosat/core/Config.cpp
    src/monosat/core/Config.h
    src/monosat/core/Dimacs.h
    src/monosat/core/Optimize.cpp
    src/monosat/core/Optimize.h
    src/monosat/core/Remap.h
    src/monosat/core/Solver.cc
    src/monosat/core/Solver.h
    src/monosat/core/SolverTypes.h
    src/monosat/core/Theory.h
    src/monosat/core/TheorySolver.h
    src/monosat/dgl/alg/dyncut/block.h
    src/monosat/dgl/alg/dyncut/graph.h
    src/monosat/dgl/alg/AugmentedSplayTree.h
    src/monosat/dgl/alg/DisjointSets.cpp
    src/monosat/dgl/alg/DisjointSets.h
    src/monosat/dgl/alg/EulerTree.h
    src/monosat/dgl/alg/Heap.h
    src/monosat/dgl/alg/LinkCut.h
    src/monosat/dgl/alg/LinkCutCost.cpp
    src/monosat/dgl/alg/LinkCutCost.h
    src/monosat/dgl/alg/SearchTree.h
    src/monosat/dgl/alg/SplayTree.h
    src/monosat/dgl/alg/Treap.h
    src/monosat/dgl/alg/TreapCustom.h
    src/monosat/dgl/AcyclicFlow.h
    src/monosat/dgl/AllPairs.cpp
    src/monosat/dgl/AllPairs.h
    src/monosat/dgl/BFS.h
    src/monosat/dgl/Chokepoint.h
    src/monosat/dgl/ConnectedComponents.cpp
    src/monosat/dgl/ConnectedComponents.h
    src/monosat/dgl/Cycle.h
    src/monosat/dgl/DFS.h
    src/monosat/dgl/DFSCycle.h
    src/monosat/dgl/Dijkstra.h
    src/monosat/dgl/DijkstraAllPairs.h
    src/monosat/dgl/Dinics.h
    src/monosat/dgl/DinicsLinkCut.h
    src/monosat/dgl/DisjointSetConnectedComponents.h
    src/monosat/dgl/Distance.cpp
    src/monosat/dgl/Distance.h
    src/monosat/dgl/DynamicConnectivity.h
    src/monosat/dgl/DynamicConnectivityImpl.h
    src/monosat/dgl/DynamicGraph.h
    src/monosat/dgl/DynamicNodes.h
    src/monosat/dgl/EdmondsKarp.h
    src/monosat/dgl/EdmondsKarpAdj.h
    src/monosat/dgl/EdmondsKarpDynamic.h
    src/monosat/dgl/FloydWarshall.h
    src/monosat/dgl/KohliTorr.h
    src/monosat/dgl/Kruskal.h
    src/monosat/dgl/MaxFlow.h
    src/monosat/dgl/MinimumSpanningTree.h
    src/monosat/dgl/NaiveDynamicConnectivity.h
    src/monosat/dgl/PKTopologicalSort.h
    src/monosat/dgl/Prim.h
    src/monosat/dgl/RamalReps.h
    src/monosat/dgl/Reach.cpp
    src/monosat/dgl/Reach.h
    src/monosat/dgl/SpiraPan.h
    src/monosat/dgl/SteinerApprox.h
    src/monosat/dgl/SteinerTree.h
    src/monosat/dgl/TarjansSCC.h
    src/monosat/dgl/ThorupDynamicConnectivity.h
    src/monosat/fsm/alg/NFAAccept.h
    src/monosat/fsm/alg/NFAGenerate.h
    src/monosat/fsm/alg/NFALinearGeneratorAcceptor.h
    src/monosat/fsm/alg/NFATransduce.h
    src/monosat/fsm/alg/NFATypes.h
    src/monosat/fsm/alg/NP0LAcceptor.h
    src/monosat/fsm/DynamicFSM.h
    src/monosat/fsm/FSMAcceptDetector.cpp
    src/monosat/fsm/FSMAcceptDetector.h
    src/monosat/fsm/FSMDetector.h
    src/monosat/fsm/FSMGeneratesDetector.cpp
    src/monosat/fsm/FSMGeneratesDetector.h
    src/monosat/fsm/FSMGeneratorAcceptorDetector.cpp
    src/monosat/fsm/FSMGeneratorAcceptorDetector.h
    src/monosat/fsm/FSMParser.h
    src/monosat/fsm/FSMTheory.h
    src/monosat/fsm/FSMTransducesDetector.cpp
    src/monosat/fsm/FSMTransducesDetector.h
    src/monosat/fsm/LSystem.h
    src/monosat/fsm/LSystemParser.h
    src/monosat/fsm/LSystemTheory.h
    src/monosat/fsm/P0LAcceptDetector.cpp
    src/monosat/fsm/P0LAcceptDetector.h
    src/monosat/geometry/alg/Heap.h
    src/monosat/geometry/bounds/BoundingBox.h
    src/monosat/geometry/bounds/BoundingSphere.h
    src/monosat/geometry/bounds/BoundingVolume.h
    src/monosat/geometry/cevans/halfspaceContainer.h
    src/monosat/geometry/cevans/halfspaceD2.h
    src/monosat/geometry/cevans/halfspaceD3.h
    src/monosat/geometry/cevans/indextable.h
    src/monosat/geometry/cevans/mathlib.cpp
    src/monosat/geometry/cevans/mathlib.h
    src/monosat/geometry/cevans/partitionspace.h
    src/monosat/geometry/cevans/point.h
    src/monosat/geometry/cevans/quickhull2D.h
    src/monosat/geometry/cevans/quickhull3D.h
    src/monosat/geometry/cevans/typedefs.h
    src/monosat/geometry/cevans/typeop.h
    src/monosat/geometry/cevans/zero.h
    src/monosat/geometry/polypartition/polypartition.cpp
    src/monosat/geometry/polypartition/polypartition.h
    src/monosat/geometry/ConvexHull.h
    src/monosat/geometry/ConvexHullCollisionDetection.cpp
    src/monosat/geometry/ConvexHullCollisionDetector.h
    src/monosat/geometry/ConvexHullDetector.cpp
    src/monosat/geometry/ConvexHullDetector.h
    src/monosat/geometry/ConvexPolygon.cpp
    src/monosat/geometry/ConvexPolygon.h
    src/monosat/geometry/Delaunay.h
    src/monosat/geometry/DelaunayPolypartition.cpp
    src/monosat/geometry/DelaunayPolypartition.h
    src/monosat/geometry/GeometryDetector.h
    src/monosat/geometry/GeometryParser.h
    src/monosat/geometry/GeometryTheory.h
    src/monosat/geometry/GeometryTypes.h
    src/monosat/geometry/GridHeightmap.cpp
    src/monosat/geometry/GridHeightmap.h
    src/monosat/geometry/Heightmap.h
    src/monosat/geometry/HeightmapDetector.h
    src/monosat/geometry/Line.h
    src/monosat/geometry/LineSegment.h
    src/monosat/geometry/MonotoneConvexHull.h
    src/monosat/geometry/PointSet.cpp
    src/monosat/geometry/PointSet.h
    src/monosat/geometry/Polygon.h
    src/monosat/geometry/PolygonSet.h
    src/monosat/geometry/QuickConvexHull.cpp
    src/monosat/geometry/QuickConvexHull.h
    src/monosat/geometry/Shape.h
    src/monosat/graph/AllPairsDetector.cpp
    src/monosat/graph/AllPairsDetector.h
    src/monosat/graph/ConnectedComponentsDetector.cpp
    src/monosat/graph/ConnectedComponentsDetector.h
    src/monosat/graph/CycleDetector.cpp
    src/monosat/graph/CycleDetector.h
    src/monosat/graph/Detector.h
    src/monosat/graph/DistanceDetector.cpp
    src/monosat/graph/DistanceDetector.h
    src/monosat/graph/GraphParser.h
    src/monosat/graph/GraphTheory.h
    src/monosat/graph/GraphTheoryTypes.h
    src/monosat/graph/MaxflowDetector.cpp
    src/monosat/graph/MaxflowDetector.h
    src/monosat/graph/MSTDetector.cpp
    src/monosat/graph/MSTDetector.h
    src/monosat/graph/ReachDetector.cpp
    src/monosat/graph/ReachDetector.h
    src/monosat/graph/SteinerDetector.cpp
    src/monosat/graph/SteinerDetector.h
    src/monosat/graph/WeightedDijkstra.h
    src/monosat/graph/WeightedDistanceDetector.cpp
    src/monosat/graph/WeightedDistanceDetector.h
    src/monosat/mtl/Alg.h
    src/monosat/mtl/Alloc.h
    src/monosat/mtl/Bitset.h
    src/monosat/mtl/Deque.h
    src/monosat/mtl/Heap.h
    src/monosat/mtl/IntMap.h
    src/monosat/mtl/IntTypes.h
    src/monosat/mtl/Map.h
    src/monosat/mtl/Queue.h
    src/monosat/mtl/Rnd.h
    src/monosat/mtl/Sort.h
    src/monosat/mtl/Vec.h
    src/monosat/mtl/XAlloc.h
    src/monosat/pb/PbParser.h
    src/monosat/pb/PbTheory.h
    src/monosat/simp/SimpSolver.cc
    src/monosat/simp/SimpSolver.h
    src/monosat/utils/Options.cc
    src/monosat/utils/Options.h
    src/monosat/utils/ParseUtils.h
    src/monosat/utils/System.cc
    src/monosat/utils/System.h
    src/monosat/Main.cc)


set(SOURCE_FILES_PB
        src/monosat/pb/Config_pb.h
        src/monosat/pb/Config_pb.cpp
        src/monosat/pb/Debug.h
        src/monosat/pb/Debug.cc
        src/monosat/pb/Hardware.h
        src/monosat/pb/Hardware_adders.cc
        src/monosat/pb/Hardware_clausify.cc
        src/monosat/pb/Hardware_sorters.cc
        src/monosat/pb/PbSolver.cc
        src/monosat/pb/PbSolver.h
        src/monosat/pb/PbSolver_convert.cc
        src/monosat/pb/PbSolver_convertAdd.cc
        src/monosat/pb/PbSolver_convertBdd.cc
        src/monosat/pb/PbSolver_convertSort.cc
        src/monosat/pb/ADTs/FEnv.h
        src/monosat/pb/ADTs/FEnv.cc
        src/monosat/pb/ADTs/File.h
        src/monosat/pb/ADTs/File.cc
        src/monosat/pb/ADTs/Global.h
        src/monosat/pb/ADTs/Global.cc
        src/monosat/pb/ADTs/Hash_standard.h
        src/monosat/pb/ADTs/Heap.h
        src/monosat/pb/ADTs/Int.h
        src/monosat/pb/ADTs/Map.h
        src/monosat/pb/ADTs/Sort.h
        src/monosat/pb/ADTs/StackAlloc.h
        src/monosat/pb/ADTs/VecAlloc.h
        src/monosat/pb/ADTs/VecMaps.h
        )

if (${PB_SUPPORT})
    set(SOURCE_FILES ${SOURCE_FILES} ${SOURCE_FILES_PB})
endif (${PB_SUPPORT})

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

add_library(libmonosat SHARED ${SOURCE_FILES} src/monosat/pb/Pb.h)
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
install(DIRECTORY src/monosat DESTINATION include
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
        PATTERN ClangShared EXCLUDE
        PATTERN GPerf EXCLUDE
        PATTERN OSX_ClangRelease EXCLUDE
        PATTERN OSX_SharedLibrary EXCLUDE
        PATTERN SharedLibrary EXCLUDE
        PATTERN SharedLibraryDbg EXCLUDE
        PATTERN OSX_SharedLibrary EXCLUDE
        PATTERN Win32 EXCLUDE
        PATTERN Win64 EXCLUDE
        PATTERN Win64SharedLibrary EXCLUDE)
