/*************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2018, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#include "monosat_MonosatJNI.h"
#include "monosat/api/Monosat.h"
#include "monosat/api/CircuitC.h"
#include <stdexcept>
#include "monosat/api/java/JNIExcept.h"

using namespace Monosat;
using namespace std;

//This code makes use of the C++/JNI exception catching
//pattern described here:
//https://stackoverflow.com/a/12014833
//Every function is wrapped in an external try/catch
//to convert any C++ exceptions into Java exceptions

JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getVersion
        (JNIEnv* env, jclass monosat_class) try{
    return env->NewStringUTF(getVersion());
}catch(...){
    javaThrow(env);
    return nullptr;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_varToLit
        (JNIEnv* env, jclass monosat_class, jint var, jboolean sign) try{
    return varToLit(var, sign);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_litToVar
        (JNIEnv* env, jclass monosat_class, jint lit) try{
    return litToVar(lit);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newSolver__
        (JNIEnv* env, jclass monosat_class) try{
    return reinterpret_cast<jlong>(newSolver());
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newSolver__Ljava_lang_String_2
        (JNIEnv* env, jclass monosat_class, jstring args) try{
    const char* str = env->GetStringUTFChars(args, 0);
    auto* ptr = newSolver_arg(str);
    env->ReleaseStringUTFChars(args, str);
    return reinterpret_cast<jlong>(ptr);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_deleteSolver
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    deleteSolver(reinterpret_cast<SolverPtr>(solverPtr));
}catch(...){
    javaThrow(env);
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_ok
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    return jboolean(ok(reinterpret_cast<SolverPtr>(solverPtr)));
}catch(...){
    javaThrow(env);
    return false;
}

JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setOutputFile
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring filename) try{
    const char* str = env->GetStringUTFChars(filename, 0);
    setOutputFile(reinterpret_cast<SolverPtr>(solverPtr), str);
    env->ReleaseStringUTFChars(filename, str);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_readGNF
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring filename) try{
    const char* str = env->GetStringUTFChars(filename, 0);
    readGNF(reinterpret_cast<SolverPtr>(solverPtr), str);
    env->ReleaseStringUTFChars(filename, str);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_loadGNF
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring filename) try{
    const char* str = env->GetStringUTFChars(filename, 0);
    loadGNF(reinterpret_cast<SolverPtr>(solverPtr), str);
    env->ReleaseStringUTFChars(filename, str);
}catch(...){
    javaThrow(env);
}

JNIEXPORT void JNICALL Java_monosat_MonosatJNI_flushFile
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    flushFile(reinterpret_cast<SolverPtr>(solverPtr));
}catch(...){
    javaThrow(env);
}

JNIEXPORT void JNICALL Java_monosat_MonosatJNI_closeFile
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    closeFile(reinterpret_cast<SolverPtr>(solverPtr));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_addLiteralName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint literal, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    addLiteralName(solver, literal, str);
    env->ReleaseStringUTFChars(name, str);
}catch(...){
    javaThrow(env);
}

JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getLiteralName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint literal, jint nameIndex) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* name = getLiteralName(solver, literal, nameIndex);
    return env->NewStringUTF(name);
}catch(...){
    javaThrow(env);
    return env->NewStringUTF("Error"); //unreachable
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_literalHasName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint literal, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    jboolean result = literalHasName(solver, literal, str);
    env->ReleaseStringUTFChars(name, str);
    return result;
}catch(...){
    javaThrow(env);
    return false; //unreachable
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_literalNameCount
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint literal) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return literalNameCount(solver, literal);
}catch(...){
    javaThrow(env);
    return false;
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_hasLiteralWithName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    bool r = hasLiteralWithName(solver, str);
    env->ReleaseStringUTFChars(name, str);
    return r;
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nNamedLiterals
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nNamedLiterals(solver));
}catch(...){
    javaThrow(env);
    return -1;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getNamedLiteralN
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint n) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    int lit = getNamedLiteralN(solver, n);
    return jint(lit);
}catch(...){
    javaThrow(env);
    return -1;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getLiteral
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    int lit = getLiteral(solver, str);
    env->ReleaseStringUTFChars(name, str);
    return jint(lit);
}catch(...){
    javaThrow(env);
    return -1;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_addVariableName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    addVariableName(solver, variable, str);
    env->ReleaseStringUTFChars(name, str);
}catch(...){
    javaThrow(env);
}

JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getVariableName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable, jint nameIndex) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* name = getVariableName(solver, variable, nameIndex);
    return env->NewStringUTF(name);
}catch(...){
    javaThrow(env);
    return env->NewStringUTF("Error"); //unreachable
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_variableHasName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    jboolean result = variableHasName(solver, variable, str);
    env->ReleaseStringUTFChars(name, str);
    return result;
}catch(...){
    javaThrow(env);
    return false; //unreachable
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_variableNameCount
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return variableNameCount(solver, variable);
}catch(...){
    javaThrow(env);
    return false;
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_hasVariableWithName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    bool r = hasVariableWithName(solver, str);
    env->ReleaseStringUTFChars(name, str);
    return r;
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nNamedVariables
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nNamedVariables(solver));
}catch(...){
    javaThrow(env);
    return -1;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nNamedBitvectors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jint(nNamedBitvectors(solver, bv));
}catch(...){
    javaThrow(env);
    return -1;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getNamedBitvectorN
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint n) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jint(getNamedBitvectorN(solver, bv, n));
}catch(...){
    javaThrow(env);
    return -1;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getNamedVariableN
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint n) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    int var = getNamedVariableN(solver, n);
    return jint(var);
}catch(...){
    javaThrow(env);
    return -1;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getVariable
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    int var = getVariable(solver, str);
    env->ReleaseStringUTFChars(name, str);
    return jint(var);
}catch(...){
    javaThrow(env);
    return -1;
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_solve
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    bool result = solve(reinterpret_cast<SolverPtr>(solverPtr));
    return jboolean(result);
}catch(...){
    javaThrow(env);
    return false;
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_solveAssumptions
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject assumptions, jint n_assumptions) try{
    return jboolean(
            solveAssumptions(reinterpret_cast<SolverPtr>(solverPtr),
                             (int*) env->GetDirectBufferAddress(assumptions),
                             n_assumptions));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setTimeLimit
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint limit) try{
    setTimeLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setConflictLimit
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint limit) try{
    setConflictLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setPropagationLimit
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint limit) try{
    setPropagationLimit(reinterpret_cast<SolverPtr>(solverPtr), limit);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_nConflicts
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    return jlong(nConflicts(reinterpret_cast<SolverPtr>(solverPtr)));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_nPropagations
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    return jlong(nPropagations(reinterpret_cast<SolverPtr>(solverPtr)));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_solveLimited
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    return solveLimited(reinterpret_cast<SolverPtr>(solverPtr));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_solveAssumptionsLimited
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject assumptions, jint n_assumptions) try{
    return solveAssumptionsLimited(reinterpret_cast<SolverPtr>(solverPtr),
                                   (int*) env->GetDirectBufferAddress(assumptions), n_assumptions);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_lastSolutionWasOptimal
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    return jboolean(lastSolutionWasOptimal(reinterpret_cast<SolverPtr>(solverPtr)));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getConflictClause
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    if(array == nullptr || length < 0){
        return getConflictClause(solver, nullptr, -1);
    }else{
        assert(length >= 0);
        return getConflictClause(solver, (int*) env->GetDirectBufferAddress(array), length);
    }
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_minimizeUnsatCore
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return minimizeUnsatCore(solver, (int*) env->GetDirectBufferAddress(array), length);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeConflictClause
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    minimizeConflictClause(solver);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_backtrack
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    backtrack(solver);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newVar
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(newVar(solver));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newNamedVar
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    jint var = jint(newNamedVar(solver, str));
    env->ReleaseStringUTFChars(name, str);
    return var;
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionVar
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable, jboolean is_decision) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionVar(solver, variable, is_decision);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_isDecisionVar
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(isDecisionVar(solver, variable));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionPriority
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable, jint priority) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionPriority(solver, variable, priority);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getDecisionPriority
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(getDecisionPriority(solver, variable));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setDecisionPolarity
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable, jboolean polarity) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    setDecisionPolarity(solver, variable, polarity);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_getDecisionPolarity
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint variable) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(getDecisionPolarity(solver, variable));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_true_1lit
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(true_lit(solver));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_disallowLiteralSimplification
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint var) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(disallowLiteralSimplification(solver, var));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_disablePreprocessing
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    disablePreprocessing(solver);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nVars
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nVars(solver));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nClauses
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nClauses(solver));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nLearnedClauses
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(nLearnedClauses(solver));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nBitvectors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jint(nBitvectors(solver, bv));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addClause
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addClause(solver, (int*) env->GetDirectBufferAddress(array), length));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addUnitClause
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addUnitClause(solver, lit));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addBinaryClause
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint l1, jint l2) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(addBinaryClause(solver, l1, l2));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_addTertiaryClause
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint l1, jint l2, jint l3) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return addTertiaryClause(solver, l1, l2, l3);
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_clearOptimizationObjectives
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    clearOptimizationObjectives(solver);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeBV
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    maximizeBV(solver, bv, bvID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeBV
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    minimizeBV(solver, bv, bvID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeLits
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    maximizeLits(solver, (int*) env->GetDirectBufferAddress(array), length);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeLits
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    minimizeLits(solver, (int*) env->GetDirectBufferAddress(array), length);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_maximizeWeightedLits
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array1, jobject array2, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    maximizeWeightedLits(solver, (int*) env->GetDirectBufferAddress(array1),
                         (int*) env->GetDirectBufferAddress(array2), length);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_minimizeWeightedLits
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array1, jobject array2, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    minimizeWeightedLits(solver, (int*) env->GetDirectBufferAddress(array1),
                         (int*) env->GetDirectBufferAddress(array2), length);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_initBVTheory
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return reinterpret_cast<jlong>(initBVTheory(solver));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector_1const
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint width, jlong constantValue) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector_const(solver, bv, width, constantValue);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector_1anon
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint width) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector_anon(solver, bv, width);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector_1lazy
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector_lazy(solver, bv, (int*) env->GetDirectBufferAddress(array), length);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBitvector
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint length) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBitvector(solver, bv, (int*) env->GetDirectBufferAddress(array), length);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_setBitvectorName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    setBitvectorName(solver, bv, bvID, str);
    env->ReleaseStringUTFChars(name, str);

}catch(...){
    javaThrow(env);
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_bitvectorHasName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return bitvectorHasName(solver, bv, bvID);
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_hasBitvectorWithName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    bool r = hasBitvectorWithName(solver, bv, str);
    env->ReleaseStringUTFChars(name, str);
    return r;
}catch(...){
    javaThrow(env);
    return false;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getBitvectorNameCount
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jint(getBitvectorNameCount(solver, bv, bvID));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getBitvectorName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, int nameIndex) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return env->NewStringUTF(getBitvectorName(solver, bv, bvID, nameIndex));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getBitvectorWidth
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return bv_width(solver, bv, bvID);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nBitvectorBits
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return bv_nBits(solver, bv, bvID);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getBitvectorBit
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jint bit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return bv_bit(solver, bv, bvID, bit);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getBitvector
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    int bvID = getBitvector(solver, bv, str);
    env->ReleaseStringUTFChars(name, str);
    return jint(bvID);
}catch(...){
    javaThrow(env);
    return -1;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1eq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_eq(solver, bv, bvID, constval);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1eq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_eq(solver, bv, bvID1, bvID2);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1neq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_neq(solver, bv, bvID, constval);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1neq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_neq(solver, bv, bvID1, bvID2);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1lt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_lt(solver, bv, bvID, constval);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1lt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_lt(solver, bv, bvID1, bvID2);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1leq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_leq(solver, bv, bvID, constval);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1leq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_leq(solver, bv, bvID1, bvID2);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1gt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_gt(solver, bv, bvID, constval);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1gt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_gt(solver, bv, bvID1, bvID2);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1const_1geq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID, jlong constval) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_const_geq(solver, bv, bvID, constval);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newBVComparison_1bv_1geq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return newBVComparison_bv_geq(solver, bv, bvID1, bvID2);
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1bitblast
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_bitblast(solver, bv, bvID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1concat
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_concat(solver, bv, aID, bID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1slice
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint lower, jint upper,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_slice(solver, bv, aID, lower, upper, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1not
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_not(solver, bv, aID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1and
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_and(solver, bv, aID, bID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1nand
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_nand(solver, bv, aID, bID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1or
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_or(solver, bv, aID, bID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1nor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_nor(solver, bv, aID, bID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1xor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_xor(solver, bv, aID, bID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1xnor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint aID, jint bID, jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_xnor(solver, bv, aID, bID, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1ite
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint condition_lit, jint bvThenID,
         jint bvElseID, jint bvResultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_ite(solver, bv, condition_lit, bvThenID, bvElseID, bvResultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1addition
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, int bvID1, int bvID2,
         int resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_addition(solver, bv, bvID1, bvID2, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1subtraction
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_subtraction(solver, bv, bvID1, bvID2, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1multiply
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_multiply(solver, bv, bvID1, bvID2, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1divide
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID1, jint bvID2,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_divide(solver, bv, bvID1, bvID2, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1min
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_min(solver, bv, (int*) env->GetDirectBufferAddress(array), n_args, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1max
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_max(solver, bv, (int*) env->GetDirectBufferAddress(array), n_args, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1popcount
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_popcount(solver, bv, (int*) env->GetDirectBufferAddress(array), n_args, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_bv_1unary
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jobject array, jint n_args,
         jint resultID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    bv_unary(solver, bv, (int*) env->GetDirectBufferAddress(array), n_args, resultID);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_at_1most_1one_1lit
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject array, jint n_args) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    at_most_one_lit(solver, (int*) env->GetDirectBufferAddress(array), n_args);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1lt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_lt(solver, rhs, n_args, (int*) env->GetDirectBufferAddress(literals),
                (int*) env->GetDirectBufferAddress(coefficients));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1leq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_leq(solver, rhs, n_args, (int*) env->GetDirectBufferAddress(literals),
                 (int*) env->GetDirectBufferAddress(coefficients));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1eq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_eq(solver, rhs, n_args, (int*) env->GetDirectBufferAddress(literals),
                (int*) env->GetDirectBufferAddress(coefficients));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1geq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_geq(solver, rhs, n_args, (int*) env->GetDirectBufferAddress(literals),
                 (int*) env->GetDirectBufferAddress(coefficients));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_assertPB_1gt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint rhs, jint n_args, jobject literals,
         jobject coefficients) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    assertPB_gt(solver, rhs, n_args, (int*) env->GetDirectBufferAddress(literals),
                (int*) env->GetDirectBufferAddress(coefficients));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_flushPB
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    flushPB(solver);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newGraph
        (JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return reinterpret_cast<jlong>(newGraph(solver));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_newGraph_1Named
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring name, jint bitwidth) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    jlong g = reinterpret_cast<jlong>(newGraph_Named(solver, str, bitwidth));
    env->ReleaseStringUTFChars(name, str);
    return g;
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getGraph
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    jlong g = reinterpret_cast<jlong>(getGraph(solver, str));
    env->ReleaseStringUTFChars(name, str);
    return g;
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getGraphWidth
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return getGraphWidth(solver, graph);
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getGraphName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return env->NewStringUTF(getGraphName(solver, graph));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jstring JNICALL Java_monosat_MonosatJNI_getNodeName
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, int nodeID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return env->NewStringUTF(getNodeName(solver, graph, nodeID));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newNode
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newNode(solver, graph));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newNode_1Named
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    jint n = jint(newNode_Named(solver, graph, str));
    env->ReleaseStringUTFChars(name, str);
    return n;
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_hasNamedNode
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jstring name) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    const char* str = env->GetStringUTFChars(name, 0);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    jboolean r = jint(hasNamedNode(solver, graph, str));
    env->ReleaseStringUTFChars(name, str);
    return r;
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newEdge
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newEdge(solver, graph, from, to, weight));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_newEdge_1bv
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(newEdge_bv(solver, graph, from, to, bvID));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nNodes
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(nNodes(solver, graph));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_nEdges
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(nEdges(solver, graph));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getEdgeLiteralN
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint n) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getEdgeLiteralN(solver, graph, n));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getEdge_1to
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint edgeLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getEdge_to(solver, graph, edgeLit));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getEdge_1from
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint edgeLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getEdge_from(solver, graph, edgeLit));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getEdge_1weight_1const
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint edgeLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getEdge_weight_const(solver, graph, edgeLit));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getEdge_1weight_1bv
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint edgeLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getEdge_weight_bv(solver, graph, edgeLit));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_edgeHasBVWeight
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint edgeLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jboolean(edgeHasBVWeight(solver, graph, edgeLit));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_reaches
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(reaches(solver, graph, from, to));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_reachesBackward
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(reachesBackward(solver, graph, from, to));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_onPath
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint nodeOnPath, jint from, jint to) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(onPath(solver, graph, nodeOnPath, from, to));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPathUnweighted_1lt_1const
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint steps) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPathUnweighted_lt_const(solver, graph, from, to, steps));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPathUnweighted_1leq_1const
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint steps) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPathUnweighted_leq_const(solver, graph, from, to, steps));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1lt_1const
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong dist) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_lt_const(solver, graph, from, to, dist));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1leq_1const
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong dist) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_leq_const(solver, graph, from, to, dist));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1lt_1bv
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_lt_bv(solver, graph, from, to, bvID));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_shortestPath_1leq_1bv
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(shortestPath_leq_bv(solver, graph, from, to, bvID));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1geq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_geq(solver, graph, from, to, weight));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1gt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jlong weight) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_gt(solver, graph, from, to, weight));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1geq_1bv
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_geq_bv(solver, graph, from, to, bvID));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_maximumFlow_1gt_1bv
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint from, jint to, jint bvID) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(maximumFlow_gt_bv(solver, graph, from, to, bvID));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_minimumSpanningTree_1leq
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(minimumSpanningTree_leq(solver, graph, weight));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_minimumSpanningTree_1lt
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(minimumSpanningTree_lt(solver, graph, weight));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_acyclic_1undirected
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(acyclic_undirected(solver, graph));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_acyclic_1directed
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(acyclic_directed(solver, graph));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_newEdgeSet
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jobject edges, jint n_edges,
         jboolean enforceEdgeAssignment) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    newEdgeSet(solver, graph, (int*) env->GetDirectBufferAddress(edges), n_edges, enforceEdgeAssignment);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_graph_1setAssignEdgesToWeight
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong weight) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    graph_setAssignEdgesToWeight(solver, graph, weight);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_createFlowRouting
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint sourceNode, jint destNode,
         jint maxflowLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return reinterpret_cast<jlong>(createFlowRouting(solver, graph, sourceNode, destNode, maxflowLit));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_addRoutingNet
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jlong routerPtr, jint disabledEdge,
         jint n_members, jobject edge_lits, jobject reach_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    FlowRouterPtr router = reinterpret_cast<FlowRouterPtr>(routerPtr);
    addRoutingNet(solver, graph, router, disabledEdge, n_members, (int*) env->GetDirectBufferAddress(edge_lits),
                  (int*) env->GetDirectBufferAddress(reach_lits));
}catch(...){
    javaThrow(env);
}

JNIEXPORT jboolean JNICALL Java_monosat_MonosatJNI_hasModel(JNIEnv* env, jclass monosat_class, jlong solverPtr) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jboolean(hasModel(solver));
}catch(...){
    javaThrow(env);
    return false;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Literal
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint literal) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(getModel_Literal(solver, literal));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getConstantModel_1Literal
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint literal) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(getConstantModel_Literal(solver, literal));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1BV
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong bitvectorPtr, jint bvID,
         jboolean getMaximumValue) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    BVTheoryPtr bv = reinterpret_cast<BVTheoryPtr>(bitvectorPtr);
    return jlong(getModel_BV(solver, bv, bvID, getMaximumValue));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1MaxFlow
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_MaxFlow(solver, graph, maxflow_literal));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1EdgeFlow
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal, jint edgeLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_EdgeFlow(solver, graph, maxflow_literal, edgeLit));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1AcyclicEdgeFlow
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint maxflow_literal, jint edgeLit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_AcyclicEdgeFlow(solver, graph, maxflow_literal, edgeLit));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jlong JNICALL Java_monosat_MonosatJNI_getModel_1MinimumSpanningTreeWeight
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint spanning_tree_literal) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jlong(getModel_MinimumSpanningTreeWeight(solver, graph, spanning_tree_literal));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1Nodes_1Length
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_Nodes_Length(solver, graph, reach_or_distance_literal));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1Nodes
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal,
         jint store_length, jobject store) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_Nodes(solver, graph, reach_or_distance_literal, store_length,
                                    (int*) env->GetDirectBufferAddress(store)));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1EdgeLits_1Length
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_EdgeLits_Length(solver, graph, reach_or_distance_literal));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_getModel_1Path_1EdgeLits
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jlong graphPtr, jint reach_or_distance_literal,
         jint store_length, jobject store) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    GraphTheorySolver_long graph = reinterpret_cast<GraphTheorySolver_long>(graphPtr);
    return jint(getModel_Path_EdgeLits(solver, graph, reach_or_distance_literal, store_length,
                                       (int*) env->GetDirectBufferAddress(store)));
}catch(...){
    javaThrow(env);
    return 0;
}

//Circuit interface



JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_And_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(And_(solver, lit_a, lit_b, lit_out));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ands_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ands_(solver, (int*) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImpliesAnd_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint implies, jobject lits, jint n_lits, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImpliesAnd_(solver, implies, (int*) env->GetDirectBufferAddress(lits), n_lits, lit_out);
}catch(...){
    javaThrow(env);
}

JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImpliesAnd
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint implies, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImpliesAnd(solver, implies, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ands
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ands(solver, (int*) env->GetDirectBufferAddress(lits), n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_And
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(And(solver, lit_a, lit_b));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Or_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Or_(solver, lit_a, lit_b, lit_out));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ors_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ors_(solver, (int*) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_ImpliesAnd
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(ImpliesAnd_(solver, (int*) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}catch(...){
    javaThrow(env);
    return 0;
}

JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_ImpliesOr
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(ImpliesOr(solver, (int*) env->GetDirectBufferAddress(lits), n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_ImpliesOr_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(ImpliesOr_(solver, (int*) env->GetDirectBufferAddress(lits), n_lits, lit_out));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImpliesOr_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint implies, jobject lits, jint n_lits, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImpliesOr_(solver, implies, (int*) env->GetDirectBufferAddress(lits), n_lits, lit_out);
}catch(...){
    javaThrow(env);
}

JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImpliesOr
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint implies, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImpliesOr(solver, implies, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ors(solver, (int*) env->GetDirectBufferAddress(lits), n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Or
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Or(solver, lit_a, lit_b));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nors(solver, (int*) env->GetDirectBufferAddress(lits), n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nor(solver, lit_a, lit_b));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nands
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nands(solver, (int*) env->GetDirectBufferAddress(lits), n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Nand
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Nand(solver, lit_a, lit_b));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xors(solver, (int*) env->GetDirectBufferAddress(lits), n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xor(solver, lit_a, lit_b));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xnors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xnors(solver, (int*) env->GetDirectBufferAddress(lits), n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Xnor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Xnor(solver, lit_a, lit_b));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Implies
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Implies(solver, lit_a, lit_b));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Implies_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Implies_(solver, lit_a, lit_b, lit_out));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ite
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_cond, jint lit_thn, jint lit_els) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ite(solver, lit_cond, lit_thn, lit_els));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Ite_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_cond, jint lit_thn, jint lit_els,
         jint lit_result) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Ite_(solver, lit_cond, lit_thn, lit_els, lit_result));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Add
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Add(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                    n_lits, (int*) env->GetDirectBufferAddress(lits_out)));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Add_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out, jint carry_lit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(Add_(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                     n_lits, (int*) env->GetDirectBufferAddress(lits_out), carry_lit));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Subtract
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(
            Subtract(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                     n_lits, (int*) env->GetDirectBufferAddress(lits_out)));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Subtract_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits,
         jobject lits_out, jint carry_lit) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(
            Subtract_(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                      n_lits, (int*) env->GetDirectBufferAddress(lits_out), carry_lit));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_Negate
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jobject lits_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    Negate(solver, (int*) env->GetDirectBufferAddress(lits), n_lits, (int*) env->GetDirectBufferAddress(lits_out));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_Negate_1
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits, jobject lits_out) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    Negate_(solver, (int*) env->GetDirectBufferAddress(lits), n_lits, (int*) env->GetDirectBufferAddress(lits_out));
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_Assert
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    Assert(solver, lit_a);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertOrTertiary
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b, jint lit_c) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertOrTertiary(solver, lit_a, lit_b, lit_c);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertOrs
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertOrs(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertOr
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertOr(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNands
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNands(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNand
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNand(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAnds
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAnds(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAnd
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAnd(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNors(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertNor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertNor(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXor(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXors(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXnors
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXnors(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertXnor
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertXnor(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertImplies
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertImplies(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertEqual
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jint lit_a, jint lit_b) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertEqual(solver, lit_a, lit_b);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAllSame
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAllSame(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_Equals
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(LEQ(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                    n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_LEQ
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(LEQ(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                    n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT jint JNICALL Java_monosat_MonosatJNI_LT
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    return jint(LT(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                   n_lits));
}catch(...){
    javaThrow(env);
    return 0;
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertEquals
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertEquals(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b),
                 n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertLEQ
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertLEQ(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertLT
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits_a, jobject lits_b, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertLT(solver, (int*) env->GetDirectBufferAddress(lits_a), (int*) env->GetDirectBufferAddress(lits_b), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertAMO
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertAMO(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}


JNIEXPORT void JNICALL Java_monosat_MonosatJNI_AssertExactlyOne
        (JNIEnv* env, jclass monosat_class, jlong solverPtr, jobject lits, jint n_lits) try{
    SolverPtr solver = reinterpret_cast<SolverPtr>(solverPtr);
    AssertExactlyOne(solver, (int*) env->GetDirectBufferAddress(lits), n_lits);
}catch(...){
    javaThrow(env);
}
