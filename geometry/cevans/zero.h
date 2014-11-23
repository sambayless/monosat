#ifndef ZERO_H
#define ZERO_H

#include "typedefs.h"
namespace cevans {
/*!
 \brief A small number.

 Numerical testing is subject to rounding errors.
 Often algorithms need to implement comparision testing
 where the in built rounding algorithms will fail.

 Different zeros for different types can be defined.
 Do not instanciate the static member in a header, instead
 instanciate it in main.cpp file so that only one definition
 is generated.

 \par Example
 \verbatim
 template<>
 double zero<double>::val = 1E-20;
 \endverbatim
 */
template<typename T>
class zero {
public:
	
	/** Generally a small positive number. */
	static T val;

	/** Test if a number is zero. */
	static boolc test(T const w) {
		if ((w < val) == false)
			return false;
		
		if ((T) 0 < w + val)
			return true;
		
		return false;
	}
	
};

}
;

#endif

