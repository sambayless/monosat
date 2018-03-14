#include "MonosatJNI.h"
#include "monosat/api/Monosat.h"

using namespace Monosat;
using namespace std;
/*
 * Class:     MonosatJNI
 * Method:    getVersion
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_MonosatJNI_getVersion
  (JNIEnv * env, jobject){
  return env->NewStringUTF(getVersion());
}

/*
 * Class:     MonosatJNI
 * Method:    varToLit
 * Signature: (IZ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_varToLit
  (JNIEnv * env, jobject, jint, jboolean){

}

/*
 * Class:     MonosatJNI
 * Method:    litToVar
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_litToVar
  (JNIEnv * env, jobject, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    newSolver
 * Signature: ()J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_newSolver__
  (JNIEnv * env, jobject){ }

/*
 * Class:     MonosatJNI
 * Method:    newSolver
 * Signature: (Ljava/lang/String;)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_newSolver__Ljava_lang_String_2
  (JNIEnv * env, jobject, jstring){ }

/*
 * Class:     MonosatJNI
 * Method:    newSolver
 * Signature: (Ljava/util/ArrayList;)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_newSolver__Ljava_util_ArrayList_2
  (JNIEnv * env, jobject, jobject){ }

/*
 * Class:     MonosatJNI
 * Method:    deleteSolver
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_deleteSolver
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    setOutputFile
 * Signature: (JLjava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setOutputFile
  (JNIEnv * env, jobject, jlong, jstring){ }

/*
 * Class:     MonosatJNI
 * Method:    readGNF
 * Signature: (JLjava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_readGNF
  (JNIEnv * env, jobject, jlong, jstring){ }

/*
 * Class:     MonosatJNI
 * Method:    solve
 * Signature: (J)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_solve
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    solveAssumptions
 * Signature: (J[I)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_solveAssumptions
  (JNIEnv * env, jobject, jlong, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    setTimeLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setTimeLimit
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    setMemoryLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setMemoryLimit
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    setConflictLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setConflictLimit
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    setPropagationLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setPropagationLimit
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    solveLimited
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_solveLimited
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    solveAssumptionsLimited
 * Signature: (J[I)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_solveAssumptionsLimited
  (JNIEnv * env, jobject, jlong, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    lastSolutionWasOptimal
 * Signature: (J)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_lastSolutionWasOptimal
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    getConflictClause
 * Signature: (J[II)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_getConflictClause
  (JNIEnv * env, jobject, jlong, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    backtrack
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_backtrack
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newVar
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newVar
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    setDecisionVar
 * Signature: (JIZ)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setDecisionVar
  (JNIEnv * env, jobject, jlong, jint, jboolean){ }

/*
 * Class:     MonosatJNI
 * Method:    isDecisionVar
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_isDecisionVar
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    setDecisionPriority
 * Signature: (JII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setDecisionPriority
  (JNIEnv * env, jobject, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getDecisionPriority
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_getDecisionPriority
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    setDecisionPolarity
 * Signature: (JIZ)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_setDecisionPolarity
  (JNIEnv * env, jobject, jlong, jint, jboolean){ }

/*
 * Class:     MonosatJNI
 * Method:    getDecisionPolarity
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_getDecisionPolarity
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    true_lit
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_true_1lit
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    disallowLiteralSimplification
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_disallowLiteralSimplification
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    disablePreprocessing
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_disablePreprocessing
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    nVars
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_nVars
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    nClauses
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_nClauses
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    nBitvectors
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_nBitvectors
  (JNIEnv * env, jobject, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    addClause
 * Signature: (J[II)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_addClause
  (JNIEnv * env, jobject, jlong, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    addUnitClause
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_addUnitClause
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    addBinaryClause
 * Signature: (JII)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_addBinaryClause
  (JNIEnv * env, jobject, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    addTertiaryClause
 * Signature: (JIII)Z
 */
JNIEXPORT jboolean JNICALL Java_MonosatJNI_addTertiaryClause
  (JNIEnv * env, jobject, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    clearOptimizationObjectives
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_clearOptimizationObjectives
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    maximizeBV
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_maximizeBV
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    minimizeBV
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_minimizeBV
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    maximizeLits
 * Signature: (J[II)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_maximizeLits
  (JNIEnv * env, jobject, jlong, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    minimizeLits
 * Signature: (J[II)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_minimizeLits
  (JNIEnv * env, jobject, jlong, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    maximizeWeightedLits
 * Signature: (J[I[II)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_maximizeWeightedLits
  (JNIEnv * env, jobject, jlong, jintArray, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    minimizeWeightedLits
 * Signature: (J[I[II)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_minimizeWeightedLits
  (JNIEnv * env, jobject, jlong, jintArray, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    initBVTheory
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_initBVTheory
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newBitvector_const
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBitvector_1const
  (JNIEnv * env, jobject, jlong, jlong, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newBitvector_anon
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBitvector_1anon
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    newBitvector
 * Signature: (JJ[II)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBitvector
  (JNIEnv * env, jobject, jlong, jlong, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_width
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_bv_1width
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_lt
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1const_1lt
  (JNIEnv * env, jobject, jlong, jlong, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_lt
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1bv_1lt
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_leq
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1const_1leq
  (JNIEnv * env, jobject, jlong, jlong, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_leq
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1bv_1leq
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_gt
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1const_1gt
  (JNIEnv * env, jobject, jlong, jlong, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_gt
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1bv_1gt
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_geq
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1const_1geq
  (JNIEnv * env, jobject, jlong, jlong, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_geq
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newBVComparison_1bv_1geq
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_bitblast
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1bitblast
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_concat
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1concat
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_slice
 * Signature: (JJIIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1slice
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_not
 * Signature: (JJII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1not
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_and
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1and
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_nand
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1nand
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_or
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1or
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_nor
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1nor
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_xor
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1xor
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_xnor
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1xnor
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_ite
 * Signature: (JJIIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1ite
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_addition
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1addition
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_subtraction
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1subtraction
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_multiply
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1multiply
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_divide
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1divide
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_min
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1min
  (JNIEnv * env, jobject, jlong, jlong, jintArray, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_max
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1max
  (JNIEnv * env, jobject, jlong, jlong, jintArray, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_popcount
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1popcount
  (JNIEnv * env, jobject, jlong, jlong, jintArray, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    bv_unary
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_bv_1unary
  (JNIEnv * env, jobject, jlong, jlong, jintArray, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    at_most_one
 * Signature: (J[II)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_at_1most_1one
  (JNIEnv * env, jobject, jlong, jintArray, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    assertPB_lt
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_assertPB_1lt
  (JNIEnv * env, jobject, jlong, jint, jint, jintArray, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    assertPB_leq
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_assertPB_1leq
  (JNIEnv * env, jobject, jlong, jint, jint, jintArray, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    assertPB_eq
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_assertPB_1eq
  (JNIEnv * env, jobject, jlong, jint, jint, jintArray, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    assertPB_geq
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_assertPB_1geq
  (JNIEnv * env, jobject, jlong, jint, jint, jintArray, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    assertPB_gt
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_assertPB_1gt
  (JNIEnv * env, jobject, jlong, jint, jint, jintArray, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    flushPB
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_flushPB
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newGraph
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_newGraph
  (JNIEnv * env, jobject, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newNode
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newNode
  (JNIEnv * env, jobject, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newEdge
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newEdge
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newEdge_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_newEdge_1bv
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    nNodes
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_nNodes
  (JNIEnv * env, jobject, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    nEdges
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_nEdges
  (JNIEnv * env, jobject, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    reaches
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_reaches
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    shortestPathUnweighted_lt_const
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_shortestPathUnweighted_1lt_1const
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    shortestPathUnweighted_leq_const
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_shortestPathUnweighted_1leq_1const
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_lt_const
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_shortestPath_1lt_1const
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_leq_const
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_shortestPath_1leq_1const
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_lt_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_shortestPath_1lt_1bv
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_leq_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_shortestPath_1leq_1bv
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_geq
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_maximumFlow_1geq
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_gt
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_maximumFlow_1gt
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_geq_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_maximumFlow_1geq_1bv
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_gt_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_maximumFlow_1gt_1bv
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    minimumSpanningTree_leq
 * Signature: (JJJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_minimumSpanningTree_1leq
  (JNIEnv * env, jobject, jlong, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    minimumSpanningTree_lt
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_minimumSpanningTree_1lt
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    acyclic_undirected
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_acyclic_1undirected
  (JNIEnv * env, jobject, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    acyclic_directed
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_acyclic_1directed
  (JNIEnv * env, jobject, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    newEdgeSet
 * Signature: (JJ[IIZ)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_newEdgeSet
  (JNIEnv * env, jobject, jlong, jlong, jintArray, jint, jboolean){ }

/*
 * Class:     MonosatJNI
 * Method:    graph_setAssignEdgesToWeight
 * Signature: (JJJ)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_graph_1setAssignEdgesToWeight
  (JNIEnv * env, jobject, jlong, jlong, jlong){ }

/*
 * Class:     MonosatJNI
 * Method:    createFlowRouting
 * Signature: (JJIII)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_createFlowRouting
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    addRoutingNet
 * Signature: (JJJII[I[I)V
 */
JNIEXPORT void JNICALL Java_MonosatJNI_addRoutingNet
  (JNIEnv * env, jobject, jlong, jlong, jlong, jint, jint, jintArray, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_Literal
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_getModel_1Literal
  (JNIEnv * env, jobject, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_BV
 * Signature: (JJIZ)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_getModel_1BV
  (JNIEnv * env, jobject, jlong, jlong, jint, jboolean){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_MaxFlow
 * Signature: (JJI)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_getModel_1MaxFlow
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_EdgeFlow
 * Signature: (JJII)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_getModel_1EdgeFlow
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_AcyclicEdgeFlow
 * Signature: (JJII)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_getModel_1AcyclicEdgeFlow
  (JNIEnv * env, jobject, jlong, jlong, jint, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_MinimumSpanningTreeWeight
 * Signature: (JJI)J
 */
JNIEXPORT jlong JNICALL Java_MonosatJNI_getModel_1MinimumSpanningTreeWeight
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_Nodes_Length
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_getModel_1Path_1Nodes_1Length
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_Nodes
 * Signature: (JJII[I)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_getModel_1Path_1Nodes
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jintArray){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_EdgeLits_Length
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_getModel_1Path_1EdgeLits_1Length
  (JNIEnv * env, jobject, jlong, jlong, jint){ }

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_EdgeLits
 * Signature: (JJII[I)I
 */
JNIEXPORT jint JNICALL Java_MonosatJNI_getModel_1Path_1EdgeLits
  (JNIEnv * env, jobject, jlong, jlong, jint, jint, jintArray){ }

