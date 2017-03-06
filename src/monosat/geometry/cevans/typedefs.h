#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <cassert>
#include <string>  
using namespace std;

typedef bool const boolc;
typedef char const charc;
typedef double const doublec;
typedef float const floatc;
typedef int const intc;
typedef string const stringc;
typedef unsigned char uchar;
typedef unsigned char const ucharc;
typedef unsigned int uint;
typedef unsigned int const uintc;
typedef unsigned long int ulongint;
typedef unsigned long int const ulongintc;
typedef unsigned short int ushortint;
typedef unsigned short int const ushortintc;

/** Display the name and value. */
#define SHOW(x) #x << '=' << (x)
#define SHOW2(x) ((x), #x )
/** Look at exact pattern. */
#define SHOW3(x) #x << "=*" << (x) << "*"
/** Get rid of flushing which is associated with a 
 stream and very complex, in my opinion endl is a very bad default in the language. */
#define endln "\n"

/** The argument is evaluated in both debug and release. */
#ifdef DEBUG_GEOMETRY
#define asserteval(xarg) xarg;
#else
#define asserteval(xarg) assert(xarg);
#endif

//  Assert Interface - debug then assert, other behaviour in release code. 
//      assertreturnOS  assert or return 1(failure).
//      assertreturnfalse  release: return false (fail)
//      assertreturnfalseN negate the argument 
//      assertreturn    assert or exit function
//      assertreturnT (xarg,retobj) user return type

/** Operating systems have 0 as success, let 1 be unsuccessful. */
#ifdef DEBUG_GEOMETRY
#define assertreturnOS(xarg) \
  assert(xarg); 
#else
#define assertreturnOS(xarg) \
{ \
  bool res = xarg;\
  if (res==false) \
    return 1; \
} 
#endif

/** If in release (not debug mode) and test fails return false. */
#ifdef DEBUG_GEOMETRY 
#define assertreturnfalse(xarg) \
assert(xarg);
#else
#define assertreturnfalse(xarg) \
{\
  bool res=(xarg);\
  if (!res) return false;\
}
#endif

/** Negate the argument. */
#define assertreturnfalseN(xarg) assertreturnfalse( ! (xarg) )

/** Assert(debug) or return(release) when the assertion fails. */
#ifdef DEBUG_GEOMETRY 
#define assertreturn(xarg) \
assert(xarg);
#else
#define assertreturn(xarg) \
{\
  bool res=(xarg);\
  if (!res) return;\
}
#endif

/** On failure user object returned. */
#ifdef DEBUG_GEOMETRY 
#define assertreturnT(xarg,failedreturnobj) \
assert(xarg);
#else
#define assertreturnT(xarg,failedreturnobj) \
{\
  bool res=(xarg);\
  if (!res) return failedreturnobj;\
}
#endif

#endif

