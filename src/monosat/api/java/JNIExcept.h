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
/*
 * Framework for converting C++ exceptions into Java exceptions
 * Derived from https://stackoverflow.com/a/12014833 (MIT License)
 */
#ifndef MONOSAT_JNIEXCEPT_H
#define MONOSAT_JNIEXCEPT_H

#include <jni.h>
#include <monosat/mtl/XAlloc.h>
#include <monosat/utils/ParseUtils.h> //for parse_error

struct ThrownJavaException : std::runtime_error {
    ThrownJavaException() : std::runtime_error(""){}

    ThrownJavaException(const std::string& msg) : std::runtime_error(msg){}
};

inline void assert_no_exception(JNIEnv* env){
    if(env->ExceptionCheck() == JNI_TRUE)
        throw ThrownJavaException("assert_no_exception");
}

//used to throw a new Java exception. Use full paths like:
//"java/lang/NoSuchFieldException"
//"java/lang/NullPointerException"
//"java/security/InvalidParameterException"
struct NewJavaException : public ThrownJavaException {
    NewJavaException(JNIEnv* env, const char* type = "", const char* message = "")
            : ThrownJavaException(type + std::string(" ") + message){
        jclass newExcCls = env->FindClass(type);
        if(newExcCls != NULL)
            env->ThrowNew(newExcCls, message);
        //if it is null, a NoClassDefFoundError was already thrown
    }
};

/**
 * Catches any C++ exception, and throws a Java exception
 * @param env The JVM environment
 */
void javaThrow(JNIEnv* env){
    try{
        throw;
    }catch(const ThrownJavaException&){
        //already reported to Java, ignore
    }catch(const std::bad_alloc& rhs){
        //translate OOM C++ exception to a Java exception
        NewJavaException(env, "java/lang/OutOfMemoryError", rhs.what());
    }catch(const Monosat::OutOfMemoryException& rhs){
        //translate OOM C++ exception to a Java exception
        NewJavaException(env, "java/lang/OutOfMemoryError", "Out of memory error");
    }catch(const std::ios_base::failure& rhs){ //sample translation
        //translate IO C++ exception to a Java exception
        NewJavaException(env, "java/io/IOException", rhs.what());
    }catch(const std::invalid_argument& e){
        NewJavaException(env, "java/lang/IllegalArgumentException", e.what());
    }catch(const Monosat::parse_error& e){
        // the natural thing to do would be to translate this into a java parse exception,
        // but that is a checked exception, which would need better handling in the JNI to properly support.
        //NewJavaException(env, "java/text/ParseException", e.what());
        // so instead, we will throw an illegal state exception
        NewJavaException(env, "java/lang/IllegalStateException", e.what());
    }catch(const std::runtime_error& e){
        //translate C++ runtime exception to a Java runtime exception
        NewJavaException(env, "java/lang/RuntimeException", e.what());
    }catch(const std::exception& e){
        //translate unknown C++ exception to a Java exception
        NewJavaException(env, "java/lang/Exception", e.what());
    }catch(...){
        //translate unknown C++ exception to a Java exception
        NewJavaException(env, "java/lang/Error", "Unknown exception type");
    }
}


#endif //MONOSAT_JNIEXCEPT_H
