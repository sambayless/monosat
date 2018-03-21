#include "monosat_MonosatJNI.h"
#include "monosat/api/Monosat.h"

using namespace Monosat;
using namespace std;
/*
 * Class:     MonosatJNI
 * Method:    getVersion
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getVersion
        (JNIEnv *env, jclass monosat_class) {
    return env->NewStringUTF(getVersion());
}

/*
 * Class:     MonosatJNI
 * Method:    varToLit
 * Signature: (IZ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_varToLit
        (JNIEnv *env, jclass monosat_class, jint var, jboolean sign) {
    return varToLit(var, sign);
}

/*
 * Class:     MonosatJNI
 * Method:    litToVar
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_litToVar
        (JNIEnv *env, jclass monosat_class, jint lit) {
    return litToVar(lit);
}

/*
 * Class:     MonosatJNI
 * Method:    newSolver
 * Signature: ()J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newSolver__
        (JNIEnv *env, jclass monosat_class) {
    return reinterpret_cast<jlong>(newSolver());
}

/*
 * Class:     MonosatJNI
 * Method:    newSolver
 * Signature: (Ljava/lang/String;)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newSolver__Ljava_lang_String_2
        (JNIEnv *env, jclass monosat_class, jstring args) {
    const char *str = env->GetStringUTFChars(args, 0);
    auto *ptr = newSolver_arg(str);
    env->ReleaseStringUTFChars(args, str);
    return reinterpret_cast<jlong>(ptr);
}

/*
 * Class:     MonosatJNI
 * Method:    deleteSolver
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_deleteSolver
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    deleteSolver(reinterpret_cast<SolverPtr>(solverPtr));

}

/*
 * Class:     MonosatJNI
 * Method:    setOutputFile
 * Signature: (JLjava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setOutputFile
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jstring filename) {
    const char *str = env->GetStringUTFChars(filename, 0);
    setOutputFile(reinterpret_cast<SolverPtr>(solverPtr), str);
    env->ReleaseStringUTFChars(filename, str);
}

/*
 * Class:     MonosatJNI
 * Method:    readGNF
 * Signature: (JLjava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_readGNF
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jstring filename) {
    const char *str = env->GetStringUTFChars(filename, 0);
    readGNF(reinterpret_cast<SolverPtr>(solverPtr), str);
    env->ReleaseStringUTFChars(filename, str);
}

/*
 * Class:     MonosatJNI
 * Method:    solve
 * Signature: (J)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_solve
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    return jboolean(solve(reinterpret_cast<SolverPtr>(solverPtr)));
}

/*
 * Class:     MonosatJNI
 * Method:    solveAssumptions
 * Signature: (J[I)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_solveAssumptions
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jintArray assumptions) {
    const jsize length = env->GetArrayLength(assumptions);
    jint *r = env->GetIntArrayElements(assumptions, NULL);
    return jboolean(solveAssumptions(reinterpret_cast<SolverPtr>(solverPtr), r, length));
}

/*
 * Class:     MonosatJNI
 * Method:    setTimeLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setTimeLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setTimeLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}

/*
 * Class:     MonosatJNI
 * Method:    setMemoryLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setMemoryLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setMemoryLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}

/*
 * Class:     MonosatJNI
 * Method:    setConflictLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setConflictLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setConflictLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}

/*
 * Class:     MonosatJNI
 * Method:    setPropagationLimit
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setPropagationLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setPropagationLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}

/*
 * Class:     MonosatJNI
 * Method:    solveLimited
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_solveLimited
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    return solveLimited(reinterpret_cast<SolverPtr>(solverPtr));
}

/*
 * Class:     MonosatJNI
 * Method:    solveAssumptionsLimited
 * Signature: (J[I)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_solveAssumptionsLimited
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jintArray assumptions) {
    const jsize length = env->GetArrayLength(assumptions);
    jint *r = env->GetIntArrayElements(assumptions, NULL);
    return solveAssumptionsLimited(reinterpret_cast<SolverPtr>(solverPtr), r, length);
}
/*
 * Class:     MonosatJNI
 * Method:    lastSolutionWasOptimal
 * Signature: (J)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_lastSolutionWasOptimal
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    return jboolean(lastSolutionWasOptimal(reinterpret_cast<SolverPtr>(solverPtr)));
}

/*
 * Class:     MonosatJNI
 * Method:    getConflictClause
 * Signature: (J[II)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getConflictClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    if (array == nullptr || length < 0) {
        return getConflictClause(solver, nullptr, -1);
    } else {
        assert(length >= 0);
        return getConflictClause(solver, (int *) env->GetDirectBufferAddress(array), length);
    }
}

/*
 * Class:     MonosatJNI
 * Method:    backtrack
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_backtrack
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    backtrack(solver);
}

/*
 * Class:     MonosatJNI
 * Method:    newVar
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newVar
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(newVar(solver));
}

/*
 * Class:     MonosatJNI
 * Method:    setDecisionVar
 * Signature: (JIZ)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionVar
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable, jboolean is_decision) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionVar(solver, variable, is_decision);
}

/*
 * Class:     MonosatJNI
 * Method:    isDecisionVar
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_isDecisionVar
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(isDecisionVar(solver, variable));
}

/*
 * Class:     MonosatJNI
 * Method:    setDecisionPriority
 * Signature: (JII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionPriority
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable, jint priority) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionPriority(solver, variable, priority);
}

/*
 * Class:     MonosatJNI
 * Method:    getDecisionPriority
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getDecisionPriority
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(getDecisionPriority(solver, variable));
}

/*
 * Class:     MonosatJNI
 * Method:    setDecisionPolarity
 * Signature: (JIZ)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionPolarity
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable, jboolean polarity) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionPolarity(solver, variable, polarity);
}

/*
 * Class:     MonosatJNI
 * Method:    getDecisionPolarity
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_getDecisionPolarity
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(getDecisionPolarity(solver, variable));
}

/*
 * Class:     MonosatJNI
 * Method:    true_lit
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_true_1lit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(true_lit(solver));
}

/*
 * Class:     MonosatJNI
 * Method:    disallowLiteralSimplification
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_disallowLiteralSimplification
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(disallowLiteralSimplification(solver, literal));
}

/*
 * Class:     MonosatJNI
 * Method:    disablePreprocessing
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_disablePreprocessing
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    disablePreprocessing(solver);
}

/*
 * Class:     MonosatJNI
 * Method:    nVars
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nVars
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nVars(solver));
}

/*
 * Class:     MonosatJNI
 * Method:    nClauses
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nClauses
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nClauses(solver));
}

/*
 * Class:     MonosatJNI
 * Method:    nBitvectors
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nBitvectors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jint(nBitvectors(solver, bv));
}

/*
 * Class:     MonosatJNI
 * Method:    addClause
 * Signature: (J[II)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addClause(solver, (int *) env->GetDirectBufferAddress(array), length));
}

/*
 * Class:     MonosatJNI
 * Method:    addUnitClause
 * Signature: (JI)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addUnitClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addUnitClause(solver, lit));
}

/*
 * Class:     MonosatJNI
 * Method:    addBinaryClause
 * Signature: (JII)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addBinaryClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint l1, jint l2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addBinaryClause(solver, l1, l2));
}

/*
 * Class:     MonosatJNI
 * Method:    addTertiaryClause
 * Signature: (JIII)Z
 */
JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addTertiaryClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint l1, jint l2, jint l3) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return addTertiaryClause(solver, l1, l2, l3);
}

/*
 * Class:     MonosatJNI
 * Method:    clearOptimizationObjectives
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_clearOptimizationObjectives
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    clearOptimizationObjectives(solver);
}

/*
 * Class:     MonosatJNI
 * Method:    maximizeBV
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeBV
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    maximizeBV(solver, bv, bvID);
}

/*
 * Class:     MonosatJNI
 * Method:    minimizeBV
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeBV
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    minimizeBV(solver, bv, bvID);
}

/*
 * Class:     MonosatJNI
 * Method:    maximizeLits
 * Signature: (J[II)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    maximizeLits(solver, (int *) env->GetDirectBufferAddress(array), length);
}

/*
 * Class:     MonosatJNI
 * Method:    minimizeLits
 * Signature: (J[II)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    minimizeLits(solver, (int *) env->GetDirectBufferAddress(array), length);
}

/*
 * Class:     MonosatJNI
 * Method:    maximizeWeightedLits
 * Signature: (J[I[II)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeWeightedLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array1, jobject array2, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    maximizeWeightedLits(solver, (int *) env->GetDirectBufferAddress(array1),
                         (int *) env->GetDirectBufferAddress(array2), length);
}

/*
 * Class:     MonosatJNI
 * Method:    minimizeWeightedLits
 * Signature: (J[I[II)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeWeightedLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array1, jobject array2, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    minimizeWeightedLits(solver, (int *) env->GetDirectBufferAddress(array1),
                         (int *) env->GetDirectBufferAddress(array2), length);
}

/*
 * Class:     MonosatJNI
 * Method:    initBVTheory
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_initBVTheory
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return reinterpret_cast<jlong>(initBVTheory(solver));
}

/*
 * Class:     MonosatJNI
 * Method:    newBitvector_const
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint width, jlong constantValue) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector_const(solver, bv, width, constantValue);
}

/*
 * Class:     MonosatJNI
 * Method:    newBitvector_anon
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector_1anon
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint width) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector_anon(solver, bv, width);
}

/*
 * Class:     MonosatJNI
 * Method:    newBitvector
 * Signature: (JJ[II)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector(solver, bv, (int *) env->GetDirectBufferAddress(array), length);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_width
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_bv_1width
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return bv_width(solver, bv, bvID);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_lt
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_lt(solver, bv, bvID, constval);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_lt
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_lt(solver, bv, bvID1, bvID2);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_leq
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_leq(solver, bv, bvID, constval);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_leq
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_leq(solver, bv, bvID1, bvID2);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_gt
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_gt(solver, bv, bvID, constval);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_gt
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_gt(solver, bv, bvID1, bvID2);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_const_geq
 * Signature: (JJIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_geq(solver, bv, bvID, constval);
}

/*
 * Class:     MonosatJNI
 * Method:    newBVComparison_bv_geq
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_geq(solver, bv, bvID1, bvID2);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_bitblast
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1bitblast
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_bitblast(solver, bv, bvID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_concat
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1concat
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_concat(solver, bv, aID, bID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_slice
 * Signature: (JJIIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1slice
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint lower, jint upper,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_slice(solver, bv, aID, lower, upper, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_not
 * Signature: (JJII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1not
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_not(solver, bv, aID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_and
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1and
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_and(solver, bv, aID, bID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_nand
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1nand
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_nand(solver, bv, aID, bID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_or
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1or
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_or(solver, bv, aID, bID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_nor
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1nor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_nor(solver, bv, aID, bID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_xor
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1xor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_xor(solver, bv, aID, bID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_xnor
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1xnor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_xnor(solver, bv, aID, bID, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_ite
 * Signature: (JJIIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1ite
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint condition_lit, jint bvThenID,
         jint bvElseID, jint bvResultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_ite(solver, bv, condition_lit, bvThenID, bvElseID, bvResultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_addition
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1addition
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, int bvID1, int bvID2, int resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_addition(solver, bv, bvID1, bvID2, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_subtraction
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1subtraction
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_subtraction(solver, bv, bvID1, bvID2, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_multiply
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1multiply
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_multiply(solver, bv, bvID1, bvID2, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_divide
 * Signature: (JJIII)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1divide
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_divide(solver, bv, bvID1, bvID2, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_min
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1min
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_min(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_max
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1max
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_max(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_popcount
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1popcount
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_popcount(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    bv_unary
 * Signature: (JJ[III)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1unary
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_unary(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}

/*
 * Class:     MonosatJNI
 * Method:    at_most_one
 * Signature: (J[II)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_at_1most_1one
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint n_args) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    at_most_one(solver, (int *) env->GetDirectBufferAddress(array), n_args);
}


/*
 * Class:     MonosatJNI
 * Method:    assertPB_lt
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_lt(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                (int *) env->GetDirectBufferAddress(coefficients));
}

/*
 * Class:     MonosatJNI
 * Method:    assertPB_leq
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_leq(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                 (int *) env->GetDirectBufferAddress(coefficients));
}

/*
 * Class:     MonosatJNI
 * Method:    assertPB_eq
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1eq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_eq(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                (int *) env->GetDirectBufferAddress(coefficients));
}

/*
 * Class:     MonosatJNI
 * Method:    assertPB_geq
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_geq(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                 (int *) env->GetDirectBufferAddress(coefficients));
}

/*
 * Class:     MonosatJNI
 * Method:    assertPB_gt
 * Signature: (JII[I[I)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_gt(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                (int *) env->GetDirectBufferAddress(coefficients));
}

/*
 * Class:     MonosatJNI
 * Method:    flushPB
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_flushPB
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    flushPB(solver);
}

/*
 * Class:     MonosatJNI
 * Method:    newGraph
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newGraph
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return reinterpret_cast<jlong>(newGraph(solver));
}

/*
 * Class:     MonosatJNI
 * Method:    newNode
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newNode
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newNode(solver,graph));
}

/*
 * Class:     MonosatJNI
 * Method:    newEdge
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newEdge
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint  to,  jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newEdge(solver,graph,from,to,weight));
}

/*
 * Class:     MonosatJNI
 * Method:    newEdge_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newEdge_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr,jint from,jint  to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newEdge(solver,graph,from,to,bvID));
}

/*
 * Class:     MonosatJNI
 * Method:    nNodes
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nNodes
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(nNodes(solver,graph));
}

/*
 * Class:     MonosatJNI
 * Method:    nEdges
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nEdges
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(nEdges(solver,graph));
}

/*
 * Class:     MonosatJNI
 * Method:    reaches
 * Signature: (JJII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_reaches
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(reaches(solver,graph,from,to));
}

/*
 * Class:     MonosatJNI
 * Method:    shortestPathUnweighted_lt_const
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPathUnweighted_1lt_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint steps) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPathUnweighted_lt_const(solver,graph,from,to,steps));
}

/*
 * Class:     MonosatJNI
 * Method:    shortestPathUnweighted_leq_const
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPathUnweighted_1leq_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint steps) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPathUnweighted_leq_const(solver,graph,from,to,steps));
}

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_lt_const
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1lt_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong dist) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_lt_const(solver,graph,from,to,dist));
}

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_leq_const
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1leq_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong dist) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_leq_const(solver,graph,from,to,dist));
}

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_lt_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1lt_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_lt_bv(solver,graph,from,to,bvID));
}

/*
 * Class:     MonosatJNI
 * Method:    shortestPath_leq_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1leq_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_leq_bv(solver,graph,from,to,bvID));
}

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_geq
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_geq(solver,graph,from,to,weight));
}

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_gt
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_gt(solver,graph,from,to,weight));
}

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_geq_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1geq_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_geq_bv(solver,graph,from,to,bvID));
}

/*
 * Class:     MonosatJNI
 * Method:    maximumFlow_gt_bv
 * Signature: (JJIII)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1gt_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_gt_bv(solver,graph,from,to,bvID));
}

/*
 * Class:     MonosatJNI
 * Method:    minimumSpanningTree_leq
 * Signature: (JJJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_minimumSpanningTree_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(minimumSpanningTree_leq(solver,graph,weight));
}

/*
 * Class:     MonosatJNI
 * Method:    minimumSpanningTree_lt
 * Signature: (JJIIJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_minimumSpanningTree_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(minimumSpanningTree_lt(solver,graph,weight));
}

/*
 * Class:     MonosatJNI
 * Method:    acyclic_undirected
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_acyclic_1undirected
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(acyclic_undirected(solver,graph));
}

/*
 * Class:     MonosatJNI
 * Method:    acyclic_directed
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_acyclic_1directed
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(acyclic_directed(solver,graph));
}

/*
 * Class:     MonosatJNI
 * Method:    newEdgeSet
 * Signature: (JJ[IIZ)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_newEdgeSet
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jobject edges, jint n_edges, jboolean enforceEdgeAssignment) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    newEdgeSet(solver,graph,(int *) env->GetDirectBufferAddress(edges),n_edges,enforceEdgeAssignment);
}

/*
 * Class:     MonosatJNI
 * Method:    graph_setAssignEdgesToWeight
 * Signature: (JJJ)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_graph_1setAssignEdgesToWeight
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    graph_setAssignEdgesToWeight(solver,graph,weight);
}

/*
 * Class:     MonosatJNI
 * Method:    createFlowRouting
 * Signature: (JJIII)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_createFlowRouting
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint sourceNode, jint destNode, jint maxflowLit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return  reinterpret_cast<jlong>(createFlowRouting(solver,graph, sourceNode, destNode, maxflowLit));
}

/*
 * Class:     MonosatJNI
 * Method:    addRoutingNet
 * Signature: (JJJII[I[I)V
 */
JNIEXPORT void JNICALL Java_monosat_MonosatJNI_addRoutingNet
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong routerPtr, jint disabledEdge, jint n_members,jobject edge_lits, jobject reach_lits){
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
     FlowRouterPtr  router = reinterpret_cast<FlowRouterPtr>(routerPtr);
    addRoutingNet(solver,graph,router,disabledEdge,n_members,(int *) env->GetDirectBufferAddress(edge_lits),(int *) env->GetDirectBufferAddress(reach_lits));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_Literal
 * Signature: (JI)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Literal
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(getModel_Literal(solver,literal));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_BV
 * Signature: (JJIZ)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1BV
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jboolean getMaximumValue) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jlong(getModel_BV(solver,bv,bvID,getMaximumValue));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_MaxFlow
 * Signature: (JJI)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1MaxFlow
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_MaxFlow(solver,graph,maxflow_literal));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_EdgeFlow
 * Signature: (JJII)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1EdgeFlow
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal, jint edgeLit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_EdgeFlow(solver,graph,maxflow_literal,edgeLit));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_AcyclicEdgeFlow
 * Signature: (JJII)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1AcyclicEdgeFlow
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal, jint edgeLit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_AcyclicEdgeFlow(solver,graph,maxflow_literal,edgeLit));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_MinimumSpanningTreeWeight
 * Signature: (JJI)J
 */
JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1MinimumSpanningTreeWeight
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint spanning_tree_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_MinimumSpanningTreeWeight(solver,graph,spanning_tree_literal));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_Nodes_Length
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1Nodes_1Length
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_Nodes_Length(solver,graph,reach_or_distance_literal));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_Nodes
 * Signature: (JJII[I)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1Nodes
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal, jint store_length, jobject store) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_Nodes(solver,graph,reach_or_distance_literal,store_length,(int *) env->GetDirectBufferAddress(store)));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_EdgeLits_Length
 * Signature: (JJI)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1EdgeLits_1Length
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_EdgeLits_Length(solver,graph,reach_or_distance_literal));
}

/*
 * Class:     MonosatJNI
 * Method:    getModel_Path_EdgeLits
 * Signature: (JJII[I)I
 */
JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1EdgeLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal,jint store_length, jobject store) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_EdgeLits(solver,graph,reach_or_distance_literal,store_length,(int *) env->GetDirectBufferAddress(store)));
}