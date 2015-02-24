/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

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

#ifndef DISTANCE_H_
#define DISTANCE_H_
#include "Reach.h"
#include <limits>
#include <vector>
#include <cstddef>
#include <gmpxx.h>

namespace dgl {

template<typename Weight>
class Distance: public Reach {
public:
	
	struct NullStatus {
		void setReachable(int u, bool reachable) {
			
		}
		bool isReachable(int u) const {
			return false;
		}
		
		void setMininumDistance(int u, bool reachable, Weight distance) {
			
		}
	};
	static NullStatus nullStatus;
	static Weight weight_unreach;
	static Weight INF;
	virtual ~Distance() {
	}
	;

	virtual void setSource(int s)=0;
	virtual int getSource()=0;
	//virtual addSource(int s)=0;
	
	virtual void update()=0;

	virtual bool connected_unsafe(int t)=0;
	virtual bool connected_unchecked(int t)=0;
	virtual bool connected(int t)=0;
	virtual Weight & unreachable() {
		return weight_unreach;
	}
	virtual Weight & inf() {
		return INF;
	}
	virtual Weight & distance(int t)=0;
	virtual Weight & distance_unsafe(int t)=0;
	virtual int previous(int node)=0;
	virtual int incomingEdge(int node)=0;
	//The maximum distance to compute up to.
	virtual void setMaxDistance(Weight & maxDistance) {
		
	}
};
template<typename Weight>
typename Distance<Weight>::NullStatus Distance<Weight>::nullStatus;

template<typename Weight>
Weight Distance<Weight>::weight_unreach = -1;

template<typename Weight>
Weight Distance<Weight>::INF = std::numeric_limits<Weight>::max()/2;

template<>
double Distance<double>::INF;
template<>
float Distance<float>::INF;

template<>
mpq_class Distance<mpq_class>::INF;
template<>
mpz_class Distance<mpz_class>::INF;
}
;

#endif /* REACH_H_ */
