#include "monosat_MonosatJNI.h"
#include "monosat/api/Monosat.h"
#include "monosat/api/CircuitC.h"

using namespace Monosat;
using namespace std;

JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getVersion
        (JNIEnv *env, jclass monosat_class) {
    return env->NewStringUTF(getVersion());
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_varToLit
        (JNIEnv *env, jclass monosat_class, jint var, jboolean sign) {
    return varToLit(var, sign);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_litToVar
        (JNIEnv *env, jclass monosat_class, jint lit) {
    return litToVar(lit);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newSolver__
        (JNIEnv *env, jclass monosat_class) {
    return reinterpret_cast<jlong>(newSolver());
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newSolver__Ljava_lang_String_2
        (JNIEnv *env, jclass monosat_class, jstring args) {
    const char *str = env->GetStringUTFChars(args, 0);
    auto *ptr = newSolver_arg(str);
    env->ReleaseStringUTFChars(args, str);
    return reinterpret_cast<jlong>(ptr);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_deleteSolver
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    deleteSolver(reinterpret_cast<SolverPtr>(solverPtr));

}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setOutputFile
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jstring filename) {
    const char *str = env->GetStringUTFChars(filename, 0);
    setOutputFile(reinterpret_cast<SolverPtr>(solverPtr), str);
    env->ReleaseStringUTFChars(filename, str);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_readGNF
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jstring filename) {
    const char *str = env->GetStringUTFChars(filename, 0);
    readGNF(reinterpret_cast<SolverPtr>(solverPtr), str);
    env->ReleaseStringUTFChars(filename, str);
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_solve
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    bool result = solve(reinterpret_cast<SolverPtr>(solverPtr));
    return jboolean(result);
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_solveAssumptions
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject assumptions, jint n_assumptions) {
    return jboolean(
            solveAssumptions(reinterpret_cast<SolverPtr>(solverPtr), (int *) env->GetDirectBufferAddress(assumptions),
                             n_assumptions));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setTimeLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setTimeLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setMemoryLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setMemoryLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setConflictLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setConflictLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setPropagationLimit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint limit) {
    setPropagationLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_solveLimited
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    return solveLimited(reinterpret_cast<SolverPtr>(solverPtr));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_solveAssumptionsLimited
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject assumptions, jint n_assumptions) {
    return solveAssumptionsLimited(reinterpret_cast<SolverPtr>(solverPtr),
                                   (int *) env->GetDirectBufferAddress(assumptions), n_assumptions);
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_lastSolutionWasOptimal
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    return jboolean(lastSolutionWasOptimal(reinterpret_cast<SolverPtr>(solverPtr)));
}


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


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_backtrack
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    backtrack(solver);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newVar
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(newVar(solver));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionVar
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable, jboolean is_decision) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionVar(solver, variable, is_decision);
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_isDecisionVar
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(isDecisionVar(solver, variable));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionPriority
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable, jint priority) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionPriority(solver, variable, priority);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getDecisionPriority
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(getDecisionPriority(solver, variable));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionPolarity
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable, jboolean polarity) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionPolarity(solver, variable, polarity);
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_getDecisionPolarity
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint variable) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(getDecisionPolarity(solver, variable));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_true_1lit
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(true_lit(solver));
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_disallowLiteralSimplification
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint var) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(disallowLiteralSimplification(solver, var));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_disablePreprocessing
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    disablePreprocessing(solver);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nVars
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nVars(solver));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nClauses
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nClauses(solver));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nBitvectors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jint(nBitvectors(solver, bv));
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addClause(solver, (int *) env->GetDirectBufferAddress(array), length));
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addUnitClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addUnitClause(solver, lit));
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addBinaryClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint l1, jint l2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addBinaryClause(solver, l1, l2));
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addTertiaryClause
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint l1, jint l2, jint l3) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return addTertiaryClause(solver, l1, l2, l3);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_clearOptimizationObjectives
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    clearOptimizationObjectives(solver);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeBV
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    maximizeBV(solver, bv, bvID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeBV
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    minimizeBV(solver, bv, bvID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    maximizeLits(solver, (int *) env->GetDirectBufferAddress(array), length);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    minimizeLits(solver, (int *) env->GetDirectBufferAddress(array), length);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeWeightedLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array1, jobject array2, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    maximizeWeightedLits(solver, (int *) env->GetDirectBufferAddress(array1),
                         (int *) env->GetDirectBufferAddress(array2), length);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeWeightedLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array1, jobject array2, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    minimizeWeightedLits(solver, (int *) env->GetDirectBufferAddress(array1),
                         (int *) env->GetDirectBufferAddress(array2), length);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_initBVTheory
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return reinterpret_cast<jlong>(initBVTheory(solver));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint width, jlong constantValue) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector_const(solver, bv, width, constantValue);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector_1anon
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint width) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector_anon(solver, bv, width);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint length) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector(solver, bv, (int *) env->GetDirectBufferAddress(array), length);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_bv_1width
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return bv_width(solver, bv, bvID);
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1eq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_eq(solver, bv, bvID, constval);
}



JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1eq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_eq(solver, bv, bvID1, bvID2);
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1neq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_neq(solver, bv, bvID, constval);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1neq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_neq(solver, bv, bvID1, bvID2);
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_lt(solver, bv, bvID, constval);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_lt(solver, bv, bvID1, bvID2);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_leq(solver, bv, bvID, constval);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_leq(solver, bv, bvID1, bvID2);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_gt(solver, bv, bvID, constval);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_gt(solver, bv, bvID1, bvID2);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_geq(solver, bv, bvID, constval);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_geq(solver, bv, bvID1, bvID2);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1bitblast
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_bitblast(solver, bv, bvID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1concat
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_concat(solver, bv, aID, bID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1slice
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint lower, jint upper,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_slice(solver, bv, aID, lower, upper, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1not
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_not(solver, bv, aID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1and
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_and(solver, bv, aID, bID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1nand
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_nand(solver, bv, aID, bID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1or
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_or(solver, bv, aID, bID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1nor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_nor(solver, bv, aID, bID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1xor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_xor(solver, bv, aID, bID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1xnor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_xnor(solver, bv, aID, bID, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1ite
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint condition_lit, jint bvThenID,
         jint bvElseID, jint bvResultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_ite(solver, bv, condition_lit, bvThenID, bvElseID, bvResultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1addition
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, int bvID1, int bvID2, int resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_addition(solver, bv, bvID1, bvID2, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1subtraction
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_subtraction(solver, bv, bvID1, bvID2, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1multiply
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_multiply(solver, bv, bvID1, bvID2, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1divide
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_divide(solver, bv, bvID1, bvID2, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1min
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_min(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1max
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_max(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1popcount
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_popcount(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1unary
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_unary(solver, bv, (int *) env->GetDirectBufferAddress(array), n_args, resultID);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_at_1most_1one
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject array, jint n_args) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    at_most_one(solver, (int *) env->GetDirectBufferAddress(array), n_args);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_lt(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                (int *) env->GetDirectBufferAddress(coefficients));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_leq(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                 (int *) env->GetDirectBufferAddress(coefficients));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1eq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_eq(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                (int *) env->GetDirectBufferAddress(coefficients));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_geq(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                 (int *) env->GetDirectBufferAddress(coefficients));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_gt(solver, rhs, n_args, (int *) env->GetDirectBufferAddress(literals),
                (int *) env->GetDirectBufferAddress(coefficients));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_flushPB
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    flushPB(solver);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newGraph
        (JNIEnv *env, jclass monosat_class, jlong solverPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return reinterpret_cast<jlong>(newGraph(solver));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newNode
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newNode(solver, graph));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newEdge
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newEdge(solver, graph, from, to, weight));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newEdge_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newEdge(solver, graph, from, to, bvID));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nNodes
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(nNodes(solver, graph));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nEdges
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(nEdges(solver, graph));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_reaches
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(reaches(solver, graph, from, to));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPathUnweighted_1lt_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint steps) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPathUnweighted_lt_const(solver, graph, from, to, steps));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPathUnweighted_1leq_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint steps) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPathUnweighted_leq_const(solver, graph, from, to, steps));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1lt_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong dist) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_lt_const(solver, graph, from, to, dist));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1leq_1const
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong dist) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_leq_const(solver, graph, from, to, dist));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1lt_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_lt_bv(solver, graph, from, to, bvID));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1leq_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_leq_bv(solver, graph, from, to, bvID));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1geq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_geq(solver, graph, from, to, weight));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1gt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_gt(solver, graph, from, to, weight));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1geq_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_geq_bv(solver, graph, from, to, bvID));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1gt_1bv
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_gt_bv(solver, graph, from, to, bvID));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_minimumSpanningTree_1leq
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(minimumSpanningTree_leq(solver, graph, weight));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_minimumSpanningTree_1lt
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(minimumSpanningTree_lt(solver, graph, weight));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_acyclic_1undirected
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(acyclic_undirected(solver, graph));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_acyclic_1directed
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(acyclic_directed(solver, graph));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_newEdgeSet
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jobject edges, jint n_edges,
         jboolean enforceEdgeAssignment) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    newEdgeSet(solver, graph, (int *) env->GetDirectBufferAddress(edges), n_edges, enforceEdgeAssignment);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_graph_1setAssignEdgesToWeight
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    graph_setAssignEdgesToWeight(solver, graph, weight);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_createFlowRouting
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint sourceNode, jint destNode,
         jint maxflowLit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return reinterpret_cast<jlong>(createFlowRouting(solver, graph, sourceNode, destNode, maxflowLit));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_addRoutingNet
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong routerPtr, jint disabledEdge,
         jint n_members, jobject edge_lits, jobject reach_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    FlowRouterPtr router = reinterpret_cast<FlowRouterPtr>(routerPtr);
    addRoutingNet(solver, graph, router, disabledEdge, n_members, (int *) env->GetDirectBufferAddress(edge_lits),
                  (int *) env->GetDirectBufferAddress(reach_lits));
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_hasModel(JNIEnv *env, jclass monosat_class, jlong solverPtr){
SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(hasModel(solver));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Literal
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(getModel_Literal(solver, literal));
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1BV
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jboolean getMaximumValue) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jlong(getModel_BV(solver, bv, bvID, getMaximumValue));
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1MaxFlow
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_MaxFlow(solver, graph, maxflow_literal));
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1EdgeFlow
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal, jint edgeLit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_EdgeFlow(solver, graph, maxflow_literal, edgeLit));
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1AcyclicEdgeFlow
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal, jint edgeLit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_AcyclicEdgeFlow(solver, graph, maxflow_literal, edgeLit));
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1MinimumSpanningTreeWeight
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint spanning_tree_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_MinimumSpanningTreeWeight(solver, graph, spanning_tree_literal));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1Nodes_1Length
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_Nodes_Length(solver, graph, reach_or_distance_literal));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1Nodes
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal,
         jint store_length, jobject store) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_Nodes(solver, graph, reach_or_distance_literal, store_length,
                                    (int *) env->GetDirectBufferAddress(store)));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1EdgeLits_1Length
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_EdgeLits_Length(solver, graph, reach_or_distance_literal));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1EdgeLits
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal,
         jint store_length, jobject store) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_EdgeLits(solver, graph, reach_or_distance_literal, store_length,
                                       (int *) env->GetDirectBufferAddress(store)));
}

//Circuit interface



JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_And_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(And_(solver, lit_a, lit_b, lit_out));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ands_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ands_(solver, (int *) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImpliesAnd_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint implies, jobject lits, jint n_lits, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImpliesAnd_(solver, implies, (int *) env->GetDirectBufferAddress(lits), n_lits, lit_out);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ands
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ands(solver, (int *) env->GetDirectBufferAddress(lits), n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_And
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(And(solver, lit_a, lit_b));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Or_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Or_(solver, lit_a, lit_b, lit_out));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ors_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ors_(solver, (int *) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_ImpliesAnd
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(ImpliesAnd(solver, (int *) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_ImpliesOr
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(ImpliesOr(solver, (int *) env->GetDirectBufferAddress(lits), n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_ImpliesOr_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(ImpliesOr_(solver, (int *) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImpliesOr_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint implies, jobject lits, jint n_lits, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImpliesOr_(solver, implies, (int *) env->GetDirectBufferAddress(lits), n_lits, lit_out);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ors(solver, (int *) env->GetDirectBufferAddress(lits), n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Or
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Or(solver, lit_a, lit_b));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nors(solver, (int *) env->GetDirectBufferAddress(lits), n_lits));

}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nor(solver, lit_a, lit_b));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nands
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nands(solver, (int *) env->GetDirectBufferAddress(lits), n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nand
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nand(solver, lit_a, lit_b));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xors(solver, (int *) env->GetDirectBufferAddress(lits), n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xor(solver, lit_a, lit_b));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xnors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xnors(solver, (int *) env->GetDirectBufferAddress(lits), n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xnor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xnor(solver, lit_a, lit_b));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Implies
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Implies(solver, lit_a, lit_b));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Implies_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Implies_(solver, lit_a, lit_b, lit_out));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ite
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_cond, jint lit_thn, jint lit_els) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ite(solver, lit_cond, lit_thn, lit_els));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ite_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_cond, jint lit_thn, jint lit_els,
         jint lit_result) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ite_(solver, lit_cond, lit_thn, lit_els, lit_result));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Add
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Add(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                    n_lits, (int *) env->GetDirectBufferAddress(lits_out)));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Add_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out, jint carry_lit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Add_(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                     n_lits, (int *) env->GetDirectBufferAddress(lits_out), carry_lit));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Subtract
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(
            Subtract(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                     n_lits, (int *) env->GetDirectBufferAddress(lits_out)));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Subtract_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out, jint carry_lit) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(
            Subtract_(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                      n_lits, (int *) env->GetDirectBufferAddress(lits_out), carry_lit));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_Negate
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jobject lits_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    Negate(solver, (int *) env->GetDirectBufferAddress(lits), n_lits, (int *) env->GetDirectBufferAddress(lits_out));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_Negate_1
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jobject lits_out) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    Negate_(solver, (int *) env->GetDirectBufferAddress(lits), n_lits, (int *) env->GetDirectBufferAddress(lits_out));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_Assert
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    Assert(solver, lit_a);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertOrTertiary
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_c) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertOrTertiary(solver, lit_a, lit_b, lit_c);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertOrs
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertOrs(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertOr
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertOr(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNands
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNands(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNand
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNand(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAnds
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAnds(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAnd
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAnd(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNors(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNor(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXor(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXors(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXnors
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXnors(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXnor
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXnor(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImplies
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImplies(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertEqual
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertEqual(solver, lit_a, lit_b);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAllSame
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAllSame(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Equals
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(LEQ(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                    n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_LEQ
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(LEQ(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                    n_lits));
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_LT
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(LT(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                   n_lits));
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertEquals
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertEquals(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b),
                 n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertLEQ
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertLEQ(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertLT
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertLT(solver, (int *) env->GetDirectBufferAddress(lits_a), (int *) env->GetDirectBufferAddress(lits_b), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAMO
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAMO(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertExactlyOne
        (JNIEnv *env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) {
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertExactlyOne(solver, (int *) env->GetDirectBufferAddress(lits), n_lits);
}
