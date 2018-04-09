
/*
 * Framework for converting C++ exceptions into Java exceptions
 * Derived from https://stackoverflow.com/a/12014833 (MIT License)
 */
#ifndef MONOSAT_JNIEXCEPT_H
#define MONOSAT_JNIEXCEPT_H
#include <jni.h>
struct ThrownJavaException : std::runtime_error {
    ThrownJavaException() :std::runtime_error("") {}
    ThrownJavaException(const std::string& msg ) :std::runtime_error(msg) {}
};
inline void assert_no_exception(JNIEnv * env) {
    if (env->ExceptionCheck()==JNI_TRUE)
        throw ThrownJavaException("assert_no_exception");
}

//used to throw a new Java exception. Use full paths like:
//"java/lang/NoSuchFieldException"
//"java/lang/NullPointerException"
//"java/security/InvalidParameterException"
struct NewJavaException : public ThrownJavaException{
    NewJavaException(JNIEnv * env, const char* type="", const char* message="")
            :ThrownJavaException(type+std::string(" ")+message)
    {
        jclass newExcCls = env->FindClass(type);
        if (newExcCls != NULL)
            env->ThrowNew(newExcCls, message);
        //if it is null, a NoClassDefFoundError was already thrown
    }
};
/**
 * Catches any C++ exception, and throws a Java exception
 * @param env The JVM environment
 */
void javaThrow(JNIEnv * env) {
    try {
        throw;
    } catch(const ThrownJavaException&) {
        //already reported to Java, ignore
    } catch(const std::bad_alloc& rhs) {
        //translate OOM C++ exception to a Java exception
        NewJavaException(env, "java/lang/OutOfMemoryError", rhs.what());
    } catch(const std::ios_base::failure& rhs) { //sample translation
        //translate IO C++ exception to a Java exception
        NewJavaException(env, "java/io/IOException", rhs.what());
    } catch(const std::exception& e) {
        //translate unknown C++ exception to a Java exception
        NewJavaException(env, "java/lang/Error", e.what());
    } catch(...) {
        //translate unknown C++ exception to a Java exception
        NewJavaException(env, "java/lang/Error", "Unknown exception type");
    }
}



#endif //MONOSAT_JNIEXCEPT_H
