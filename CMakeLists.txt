cmake_minimum_required(VERSION 3.02)
project(monosat)
include(GNUInstallDirs)
#enable to include built-in support for Pseudo-Boolean constraints (recommended)
option (PB_SUPPORT "Combile with Pseudo-Boolean support (using Minisat+)" ON)
#set to compile Java JNI support into the shared library (requires a JDK)
option(JAVA "Build the Java library" OFF)
#set to install python support (requires python 3+)
option(PYTHON "Install the Python3 library" OFF)
#set to install python using the cython bindings (only takes effect if PYTHON is set to ON)
option(CYTHON "Use Cython for Python bindings (requires Cython)" OFF)
#set to OFF to disable linking GPL sources
option(GPL "Link GPLv2 sources, so that the compiled binary is licensed under the terms of the GPLv2, rather than MIT (greatly improves the performance of maximum flow predicates significantly)." ON)
option(SHOW_GIT_VERSION "Include git --describe in the build version" ON)
set (JAVA_SOURCE_FILES "")
set (JAVA_NATIVE_SOURCE_FILES "")

#set(CMAKE_ENABLE_COMPILE_COMMANDS ON) # set to enable clang tidy
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DNO_GMP -std=c++11 -Werror=return-type -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare  -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-warning-option -Werror=return-stack-address")
else()
    # compiler is gcc
endif()


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -DDEBUG_MONOSAT -DDEBUG_CORE " ) #  -DDEBUG_GRAPH -DDEBUG_BV -DDEBUG_THEORY -DDEBUG_DGL -DDEBUG_GRAPH -DDEBUG_CORE -DEBUG_GEOMETRY -DEBUG_FSM -DEBUG_PB
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3") #it is _critically_ important that release builds specify -DNDEBUG!
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O2") #it is _critically_ important that release builds specify -DNDEBUG!
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DNDEBUG") #it is _critically_ important that release builds specify -DNDEBUG!

set( CMAKE_VERBOSE_MAKEFILE ON )

#enable more flexible shared library path loading on OSX
set(CMAKE_MACOSX_RPATH 1)

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
        src/monosat/dgl/DynamicBackGraph.h
        src/monosat/dgl/Graph.h
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


if (JAVA)
    find_package(Java 1.8 REQUIRED)
    find_package(JNI 1.8 REQUIRED)
    include(UseJava)

    set(JAVA_SOURCE_FILES ${JAVA_SOURCE_FILES}
            src/monosat/api/java/monosat/MonosatJNI.java
            src/monosat/api/java/monosat/Solver.java
            src/monosat/api/java/monosat/Logic.java
            src/monosat/api/java/monosat/Lit.java
            src/monosat/api/java/monosat/Comparison.java
            src/monosat/api/java/monosat/BitVector.java
            src/monosat/api/java/monosat/Graph.java
            src/monosat/api/java/monosat/NoModelException.java
            src/monosat/api/java/monosat/TrivialContradictionException.java
            )

    message(STATUS "Compiling shared library with support for JNI")
    set_source_files_properties(src/monosat/api/java/MonosatJNI.cpp OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/monosat/api/java/monosat_MonosatJNI.h)
    #set_source_files_properties(src/monosat/api/java/monosat_MonosatJNI.h OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/__MonosatJNI.h)
    set_source_files_properties(src/monosat/api/java/monosat_MonosatJNI.h OBJECT_DEPENDS "${JAVA_SOURCE_FILES}")
    # generate MonosatJNI.h header
    # consider using 'create_javah' instead of custom command
    add_custom_command(
            PRE_BUILD
            #OUTPUT ${CMAKE_SOURCE_DIR}/src/monosat/api/java/monosat_MonosatJNI.h
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/__MonosatJNI.h # Using a fake file here that is never actually generated, to ensure that javac is always re-run
            #javah is deprecated, using javac -h instead
            # COMMAND ${Java_JAVAH_EXECUTABLE} -jni -verbose  monosat.MonosatJNI
            #COMMAND echo "Running ${Java_JAVAC_EXECUTABLE} from ${CMAKE_SOURCE_DIR}/src/monosat/api/java/"
            #COMMAND echo "${Java_JAVAC_EXECUTABLE}"
            #COMMAND ${Java_JAVAC_EXECUTABLE} -d ${CMAKE_CURRENT_BINARY_DIR}/src/monosat/api/java/monosat  monosat/MonosatJNI.java
            COMMAND ${Java_JAVAC_EXECUTABLE} -h ${CMAKE_SOURCE_DIR}/src/monosat/api/java monosat/MonosatJNI.java
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/monosat/api/java/
    )

    set(JAVA_NATIVE_SOURCE_FILES ${JAVA_NATIVE_SOURCE_FILES} src/monosat/api/java/MonosatJNI.cpp src/monosat/api/java/monosat_MonosatJNI.h)

    if (JNI_FOUND)
        message (STATUS "JNI_INCLUDE_DIRS=${JNI_INCLUDE_DIRS}")
        message (STATUS "JNI_LIBRARIES=${JNI_LIBRARIES}")
    endif()
    include_directories(${JNI_INCLUDE_DIRS})
    set(CMAKE_JNI_TARGET TRUE)



else()
    message(STATUS "Compiling wihtout library support for Java. To enable Java support, set -DJAVA=ON")
endif (JAVA)

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




add_library(libmonosat_static STATIC ${SOURCE_FILES})
set_target_properties(libmonosat_static PROPERTIES OUTPUT_NAME monosat)
target_link_libraries(libmonosat_static z.a)
#target_link_libraries(libmonosat_static m.a) # c++ doesn't require libm to be explicitly linked
target_link_libraries(libmonosat_static gmpxx.a)
target_link_libraries(libmonosat_static gmp.a)

if (UNIX)
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    if (NOT APPLE)
        target_link_libraries(libmonosat_static rt)
    else()
        target_link_libraries(libmonosat_static)
    endif()
endif (UNIX)


#add_executable(monosat_static ${SOURCE_FILES})
add_executable(monosat_static src/monosat/Main.cc)

set_target_properties(monosat_static PROPERTIES OUTPUT_NAME monosat)
if (NOT APPLE)
    target_link_libraries (monosat_static "-static-libgcc -static-libstdc++")
endif()
target_link_libraries(monosat_static libmonosat_static)
target_link_libraries(monosat_static z.a)
#target_link_libraries(monosat_static m.a)  # c++ doesn't require libm to be explicitly linked
target_link_libraries(monosat_static gmpxx.a)
target_link_libraries(monosat_static gmp.a)



if (UNIX)
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    if (NOT APPLE)
        target_link_libraries(monosat_static rt)
    else()
        target_link_libraries(monosat_static)
    endif()
endif (UNIX)



add_library(libmonosat SHARED ${SOURCE_FILES} ${JAVA_NATIVE_SOURCE_FILES})
set_target_properties(libmonosat PROPERTIES OUTPUT_NAME monosat)
target_link_libraries(libmonosat z)
if (NOT APPLE)
    #target_link_libraries(libmonosat m)  # c++ doesn't require libm to be explicitly linked
endif()
target_link_libraries(libmonosat gmpxx)
target_link_libraries(libmonosat gmp)
if (JAVA)
    target_link_libraries(libmonosat ${JNI_LIBRARIES})
endif (JAVA)
if (UNIX)
    #librt is needed for clock_gettime, which is enabled for linux only
    #(clock_gettime is used for capturing detailed timing statistics only)
    if (NOT APPLE)
        target_link_libraries(libmonosat rt)
    else()
        target_link_libraries(libmonosat)
    endif()
endif (UNIX)


add_executable(monosat  src/monosat/Main.cc)
#By default, build the staticly linked version of monosat rather than this dynamically linked one.
set_target_properties(monosat PROPERTIES EXCLUDE_FROM_ALL 1)
target_link_libraries(monosat libmonosat)
target_link_libraries(monosat z)
if (NOT APPLE)
   # target_link_libraries(monosat m) # c++ doesn't require libm to be explicitly linked
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

if (JAVA)
    target_link_libraries(libmonosat ${JNI_LIBRARIES}) #Not clear if this is required

    add_jar(monosat_jar ${JAVA_SOURCE_FILES} OUTPUT_NAME monosat)
    get_target_property(_jarFile monosat_jar JAR_FILE)
    get_target_property(_classDir monosat_jar CLASSDIR)
    message(STATUS "Monosat jar file will be created at ${CMAKE_BINARY_DIR}/${_jarFile}")


    file(GLOB_RECURSE RELATIVE_JAVA_PATHS RELATIVE "${CMAKE_SOURCE_DIR}/src/monosat/api/java/" ${JAVA_SOURCE_FILES} )

    #message(STATUS "Monosat javadoc will be compiled to ${CMAKE_BINARY_DIR}/javadoc, and ${CMAKE_BINARY_DIR}/monosat_javadoc.zip")

    message(STATUS "Monosat java source will be stored in to ${CMAKE_BINARY_DIR}/monosat_src.zip")

    #Appending  "|| true" to each of the commands below, to suppress non-zero exit codes in case of errors,
    #To ensure that these commands, if they fail, cannot prevent the libmonosat target as a whole from building.
    #Is there a better way to make custom commands non-failing/optional in cmake?
    #This will likely only work correctly for unix targets
    add_custom_command(TARGET libmonosat
            POST_BUILD
            COMMAND
            "${Java_JAR_EXECUTABLE}" "cfM" "${CMAKE_BINARY_DIR}/monosat_src.zip" ${RELATIVE_JAVA_PATHS}  "||" "true" # do not quote RELATIVE_JAVA_PATHS, or else semicolons will be injected between items.
            #OUTPUT "${CMAKE_BINARY_DIR}/monosat_src.zip"
            WORKING_DIRECTORY  "${CMAKE_SOURCE_DIR}/src/monosat/api/java/"
            )

    #add_custom_command(TARGET libmonosat
    #        POST_BUILD
    #        COMMAND
    #        "${Java_JAVADOC_EXECUTABLE}" "-d" "${CMAKE_CURRENT_BINARY_DIR}/javadoc"  ${RELATIVE_JAVA_PATHS} "||" "true" # do not quote RELATIVE_JAVA_PATHS, or else semicolons will be injected between items.
    #        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/src/monosat/api/java/"
    #        #OUTPUT "${CMAKE_BINARY_DIR}/javadoc"
    #        )
    #add_custom_command(TARGET libmonosat
    #        POST_BUILD
    #        DEPENDS ${CMAKE_BINARY_DIR}/javadoc
    #        COMMAND
    #        "${Java_JAR_EXECUTABLE}" "cfM" "monosat_javadoc.zip" "-C" "${CMAKE_CURRENT_BINARY_DIR}" "javadoc"  "||" "true"
    #        #OUTPUT "${CMAKE_BINARY_DIR}/monosat_javadoc.zip"
    #        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    #        )
    #add_custom_command(TARGET libmonosat
    #        POST_BUILD
    #        DEPENDS ${CMAKE_BINARY_DIR}/javadoc
    #        COMMAND cmake -E remove_directory  "${CMAKE_BINARY_DIR}/javadoc" "||" "true"
    #        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    #        )
    #ensure the source, doc archives are cleaned with make clean
    SET_DIRECTORY_PROPERTIES(PROPERTIES
            ADDITIONAL_MAKE_CLEAN_FILES "monosat_src.zip")

    add_dependencies(libmonosat monosat_jar)

endif (JAVA)


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

if (JAVA)

    #set(JAVA_LIB_INSTALL_DIR /usr/${CMAKE_INSTALL_LIBDIR}/java)
    message(STATUS "Installing jar to ${CMAKE_INSTALL_FULL_LIBDIR}")
    install_jar(monosat_jar ${CMAKE_INSTALL_FULL_LIBDIR})
    #install_jni_symlink(monosat_jar . )
endif(JAVA)

if (PYTHON)
    message(STATUS "Including Python library support")
    set(Python_ADDITIONAL_VERSIONS 3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8)
    # following https://bloerg.net/2012/11/10/cmake-and-distutils.html
    find_package(PythonInterp 3 REQUIRED)

    set(SETUP_PY "${CMAKE_BINARY_DIR}/setup.py")
    set(DEPS      "${CMAKE_SOURCE_DIR}/src/monosat/api/python/monosat/__init__.py")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")
    set(PACKAGE_DIR "${CMAKE_SOURCE_DIR}/src/monosat/api/python/")
    if (CYTHON)
        message(STATUS "Python library will use Cython (if available)")
        set(USE_CYTHON "True")
    endif()

    configure_file("${PACKAGE_DIR}/setup.py" ${SETUP_PY})
    add_custom_command(OUTPUT ${OUTPUT}
            COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build
            COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
            DEPENDS ${DEPS})

    add_custom_target(target ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install -f --prefix=${CMAKE_INSTALL_PREFIX})")
else()
    message(STATUS "Not including Python library support. To enable Python support, set -DPYTHON=ON.")
endif()