/**************************************************************************************************
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
#ifndef DELAUNAYPOLYPARTITION_H_
#define DELAUNAYPOLYPARTITION_H_

#include "Delaunay.h"

#include "PolygonSet.h"
#include <gmpxx.h>
#include <vector>

//interface to the polypartition library
template<unsigned int D, class T>
class DelaunayPolypartition: public Delaunay<D, T> {
	
	PolygonSet<D, T> & polygons;
	std::vector<ConvexPolygon<D, T> > triangulation;

public:
	DelaunayPolypartition(PolygonSet<D, T> & p) :
			polygons(p) {
		
	}
	
	void update();

	std::vector<ConvexPolygon<D, T> > & getTriangulation() {
		update();
		return triangulation;
	}
	
private:
	
};
template<>
DelaunayPolypartition<2, double>::DelaunayPolypartition(PolygonSet<2, double> & p);
template<>
DelaunayPolypartition<2, mpq_class>::DelaunayPolypartition(PolygonSet<2, mpq_class> & p);
template<>
void DelaunayPolypartition<2, double>::update();
template<>
void DelaunayPolypartition<2, mpq_class>::update();
#endif /* DELAUNAYPOLYPARTITION_H_ */
