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

#include "DelaunayPolypartition.h"

#include "polypartition/polypartition.h"
template<>
DelaunayPolypartition<2, double>::DelaunayPolypartition(PolygonSet<2, double> & p) :
		polygons(p) {
	
}
template<>
DelaunayPolypartition<2, mpq_class>::DelaunayPolypartition(PolygonSet<2, mpq_class> & p) :
		polygons(p) {
	
}

template<>
void DelaunayPolypartition<2, double>::update() {
	TPPLPartition<double> pp;
	
	list<TPPLPoly<double>> input_polygons, result;
	
	for (int i = 0; i < polygons.size(); i++) {
		
		if (polygons.isEnabled(i)) {
			auto & polygon = polygons[i];
			input_polygons.push_back(TPPLPoly<double>());
			TPPLPoly<double> & poly = input_polygons.back();
			poly.Init(polygon.size());
			for (auto & p : polygon) {
				poly[i].x = p.x;
				poly[i].y = p.y;
			}
		}
	}
	
	pp.Triangulate_MONO(&input_polygons, &result);
	
}

template<>
void DelaunayPolypartition<2, mpq_class>::update() {
	TPPLPartition<mpq_class> pp;
	
	list<TPPLPoly<mpq_class>> input_polygons, result;
	
	for (int i = 0; i < polygons.size(); i++) {
		
		if (polygons.isEnabled(i)) {
			auto & polygon = polygons[i];
			input_polygons.push_back(TPPLPoly<mpq_class>());
			TPPLPoly<mpq_class> & poly = input_polygons.back();
			poly.Init(polygon.size());
			for (auto & p : polygon) {
				poly[i].x = p.x;
				poly[i].y = p.y;
			}
		}
	}
	
	pp.Triangulate_MONO(&input_polygons, &result);
	
}

