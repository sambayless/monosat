cmake_minimum_required(VERSION 2.7)
project(monosat)

#enable to include built-in support for Pseudo-Boolean constraints (recommended)
option (PB_SUPPORT "Combile with Pseudo-Boolean support (using Minisat+)" ON)
#set to compile Java JNI support into the shared library (requires a JDK)
option(JAVA_SUPPORT "Build the static library" OFF)
#set to OFF to disable linking GPL sources
option(GPL "Link GPLv2 sources, so that the compiled binary is licensed under the terms of the GPLv2, rather than MIT (greatly improves the performance of maximum flow predicates significantly)." ON)
option(SHOW_GIT_VERSION "Include git --describe in the build version" ON)
set (JAVA_SOURCE_FILES "")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DNO_GMP -std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare  -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -DDEBUG_MONOSAT -DDEBUG_CORE " ) #  -DDEBUG_GRAPH -DDEBUG_BV -DDEBUG_THEORY -DDEBUG_DGL -DDEBUG_GRAPH -DDEBUG_CORE -DEBUG_GEOMETRY -DEBUG_FSM -DEBUG_PB
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3") #it is _critically_ important that release builds specify -DNDEBUG!
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O2") #it is _critically_ important that release builds specify -DNDEBUG!
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DNDEBUG") #it is _critically_ important that release builds specify -DNDEBUG!

set( CMAKE_VERBOSE_MAKEFILE on )

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.\nUse cmake -DCMAKE_BUILD_TYPE=Debug to select Debug release.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"  "MinSizeRel" "RelWithDebInfo" "Sanitize")
endif()

MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} " (set with -DCMAKE_BUILD_TYPE=<Release/Debug>)" )

message("* Adding build types...")
SET(CMAKE_CXX_FLAGS_SANITIZE
        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used by the C++ compiler during address sanitize  builds."
        FORCE )
SET(CMAKE_C_FLAGS_SANITIZE
        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used by the C compiler during address sanitize builds."
        FORCE )
SET(CMAKE_EXE_LINKER_FLAGS_SANITIZE
        "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used for linking binaries during address sanitize builds."
        FORCE )
SET(CMAKE_SHARED_LINKER_FLAGS_SANITIZE
        "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined"
        CACHE STRING "Flags used by the shared libraries linker during address sanitize builds."
        FORCE )
MARK_AS_ADVANCED(
        CMAKE_CXX_FLAGS_SANETIZE
        CMAKE_C_FLAGS_SANETIZE
        CMAKE_EXE_LINKER_FLAGS_SANETIZE
        CMAKE_SHARED_LINKER_FLAGS_SANETIZE )


include_directories(src/)

set(SOURCE_FILES
        src/monosat/Version.h
        src/monosat/Version.cpp
        src/monosat/amo/AMOParser.h
        src/monosat/amo/AMOTheory.h
        src/monosat/api/Circuit.h
        src/monosat/api/Monosat.cpp
        src/monosat/api/Monosat.h
        src/monosat/api/CircuitC.h
        src/monosat/api/CircuitC.cpp
        src/monosat/api/MonosatInternal.h
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
        src/monosat/dgl/RamalRepsBatched.h
        src/monosat/dgl/RamalRepsBatchedUnified.h
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
        src/monosat/fsm/alg/NFAAcceptor.h
        src/monosat/fsm/alg/NFAGraphAccept.h
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
        src/monosat/Main.cc
        src/monosat/core/Heuristic.h
        src/monosat/api/Logic.h
        src/monosat/graph/GraphHeuristic.h
        src/monosat/dgl/CachedReach.h
        src/monosat/routing/FlowRouter.h
        src/monosat/routing/FlowRouterParser.h
        )

#GPL Licenced source files are only included if the MIT option is not set
set(SOURCE_FILES_GPL
        src/monosat/dgl/alg/dyncut/block.h
        src/monosat/dgl/alg/dyncut/graph.h
)

set(SOURCE_FILES_PB
        src/monosat/pb/Pb.h
        src/monosat/pb/Clausify.h
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

if (PB_SUPPORT)
    set(SOURCE_FILES ${SOURCE_FILES} ${SOURCE_FILES_PB})
endif (PB_SUPPORT)


if (JAVA_SUPPORT)
    find_package(Java 1.8 REQUIRED)
    find_package(JNI 1.8 REQUIRED)

    message(STATUS "Compiling shared library with support for JNI")
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/src/monosat/api/java/MonosatJNI.cpp DEPENDS ${CMAKE_SOURCE_DIR}/src/monosat/api/java/monosat_MonosatJNI.h)

    # generate MonosatJNI.h header
    add_custom_command(
            PRE_BUILD
            OUTPUT ${CMAKE_SOURCE_DIR}/src/monosat/api/java/monosat_MonosatJNI.h
            #javah is now deprecated, using javac -h instead
            # COMMAND ${Java_JAVAH_EXECUTABLE} -jni -verbose  monosat.MonosatJNI
            #COMMAND echo "Running ${Java_JAVAC_EXECUTABLE} from ${CMAKE_SOURCE_DIR}/src/monosat/api/java/"
            #COMMAND echo "${Java_JAVAC_EXECUTABLE}"
            #COMMAND ${Java_JAVAC_EXECUTABLE} -d ${CMAKE_CURRENT_BINARY_DIR}/src/monosat/api/java/monosat  monosat/MonosatJNI.java
            COMMAND ${Java_JAVAC_EXECUTABLE} -h ${CMAKE_SOURCE_DIR}/src/monosat/api/java monosat/MonosatJNI.java
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/monosat/api/java/
    )

    set(JAVA_SOURCE_FILES ${JAVA_SOURCE_FILES} src/monosat/api/java/MonosatJNI.cpp src/monosat/api/java/monosat_MonosatJNI.h)

    if (JNI_FOUND)
        message (STATUS "JNI_INCLUDE_DIRS=${JNI_INCLUDE_DIRS}")
        message (STATUS "JNI_LIBRARIES=${JNI_LIBRARIES}")
    endif()
    include_directories(${JNI_INCLUDE_DIRS})
endif (JAVA_SUPPORT)

if (GPL)
    MESSAGE( STATUS "Linking GPLv2 source files. Use \"cmake -DGPL=OFF\" to build without GPL sources." )
    add_definitions(-DLINK_GPL)
    set(SOURCE_FILES ${SOURCE_FILES} ${SOURCE_FILES_GPL})
endif()


if (SHOW_GIT_VERSION)
    #only attempt to load the git version if the git repo is found in the source directory
    if( EXISTS "${CMAKE_SOURCE_DIR}/.git" )
        MESSAGE( STATUS "Found git repository" )
        SET(_definitions -DMONOSAT_VERSION=\"$CMAKE_BUILD_TYPE $$\(shell git describe --abbrev=4 --dirty --always --tags\)\")
        set_source_files_properties(${CMAKE_SOURCE_DIR}/src/monosat/Version.cpp VERBATIM PROPERTIES COMPILE_FLAGS "-DMONOSAT_BUILD=${CMAKE_BUILD_TYPE} -DMONOSAT_VERSION=\"$\(shell git describe --abbrev=4 --dirty --always --tags\)\"")
        set_source_files_properties(${CMAKE_SOURCE_DIR}/src/monosat/Version.cpp DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/__Version.cpp)
        add_custom_command(
                OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/__Version.cpp #fake file to force that Version.cpp is always rebuilt
                COMMAND echo
                ARGS \"MonoSAT ${CMAKE_BUILD_TYPE} build, version: $\(shell git describe --abbrev=4 --dirty --always --tags\)\"
        )
    endif()
endif (SHOW_GIT_VERSION)

add_executable(monosat ${SOURCE_FILES})
#set_target_properties(monosat PROPERTIES EXCLUDE_FROM_ALL 1)
target_link_libraries(monosat z)
if (NOT APPLE)
    target_link_libraries(monosat m)
endif()
target_link_libraries(monosat gmpxx)
target_link_libraries(monosat gmp)


if (UNIX)
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    if (NOT APPLE)
        target_link_libraries(monosat rt)
    else()
        target_link_libraries(monosat)
    endif()
endif (UNIX)

#add_executable(monosat_static ${SOURCE_FILES})
#set_target_properties(monosat_static PROPERTIES OUTPUT_NAME monosat)
#if (NOT APPLE)
#    target_link_libraries (monosat_static "-static-libgcc -static-libstdc++")
#endif()
#target_link_libraries(monosat_static z.a)
#target_link_libraries(monosat_static m.a)
#target_link_libraries(monosat_static gmpxx.a)
#target_link_libraries(monosat_static gmp.a)



#if (UNIX)
#    #librt is needed for clock_gettime, which is enabled for linux only
#    #(clock_gettime is used for capturing detailed timing statistics only)
#    if (NOT APPLE)
#        target_link_libraries(monosat_static rt)
#    else()
#        target_link_libraries(monosat_static)
#    endif()
#endif (UNIX)



add_library(libmonosat SHARED ${SOURCE_FILES} ${JAVA_SOURCE_FILES})
set_target_properties(libmonosat PROPERTIES OUTPUT_NAME monosat)
target_link_libraries(libmonosat z)
if (NOT APPLE)
    target_link_libraries(libmonosat m)
endif()
target_link_libraries(libmonosat gmpxx)
target_link_libraries(libmonosat gmp)
if (JAVA_SUPPORT)
    target_link_libraries(libmonosat ${JNI_LIBRARIES})
endif (JAVA_SUPPORT)
if (UNIX)
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    if (NOT APPLE)
        target_link_libraries(libmonosat rt)
    else()
        target_link_libraries(libmonosat)
    endif()
endif (UNIX)


#add_library(libmonosat_static ${SOURCE_FILES})
#set_target_properties(libmonosat_static PROPERTIES OUTPUT_NAME monosat)
#target_link_libraries(libmonosat_static z.a)
#target_link_libraries(libmonosat_static m.a)
#target_link_libraries(libmonosat_static gmpxx.a)
#target_link_libraries(libmonosat_static gmp.a)

#if (UNIX)
#    #librt is needed for clock_gettime, which is enabled for linux only
#    #(clock_gettime is used for capturing detailed timing statistics only)
#    if (NOT APPLE)
#        target_link_libraries(libmonosat_static rt)
#    else()
#        target_link_libraries(libmonosat_static)
#    endif()
#endif (UNIX)

install(TARGETS monosat libmonosat #libmonosat_static
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
        PATTERN Test EXCLUDE
        PATTERN TestRelease EXCLUDE
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