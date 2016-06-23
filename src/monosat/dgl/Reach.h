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

#ifndef REACH_H_
#define REACH_H_

#include <vector>
namespace dgl {

class Reach {
public:
	
	struct NullStatus {
		void setReachable(int u, bool reachable) {
			
		}
		bool isReachable(int u) const {
			return false;
		}

		void setMininumDistance(int u, bool reachable, int distance) {

		}

	};
	static NullStatus nullStatus;

	virtual int numUpdates() const=0;
	virtual ~Reach() {
	}
	;

	virtual void setSource(int s)=0;
	virtual int getSource()=0;
	//virtual addSource(int s)=0;
	
	virtual void update()=0;

	virtual bool connected_unsafe(int t)=0;
	virtual bool connected_unchecked(int t)=0;
	virtual bool connected(int t)=0;
	//virtual int distance( int t)=0;
	//virtual int distance_unsafe(int t)=0;
	virtual int previous(int node)=0;
	virtual int incomingEdge(int node)=0;
	//The maximum distance to compute up to.
	/*	virtual void setMaxDistance(int maxDistance){

	 }*/
};
}
;
#endif /* REACH_H_ */
