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
#ifndef POLYGONSET_H_
#define POLYGONSET_H_
#include <vector>
#include "GeometryTypes.h"
/**
 * A dynamic point set
 */
template<unsigned int D, class T = double>
class PolygonSet {
	//Dimension of this shape
	std::vector<Polygon<D, T>> polygons;
	std::vector<bool> enabled;

	int num_enabled = 0;

public:
	
	int dimension() {
		return D;
	}
	int size() const {
		return polygons.size();
	}
	
	int nEnabled() const {
		return num_enabled;
	}
	bool isEnabled(int polyID) {
		return enabled[polyID];
	}
	
	void setPolygonEnabled(int polyID, bool _enabled) {
		if (isEnabled(polyID) == _enabled)
			return;
		enabled[polyID] = _enabled;
		if (_enabled) {
			num_enabled++;
		} else {
			num_enabled--;
		}
	}
	
	Point<D, T> & getPoint(int polyID, int pointID) {
		return polygons[polyID][pointID];
	}
	
	std::vector<Polygon<D, T>> & getEnabledPolygons(std::vector<Polygon<D, T>> & points_out) {
		points_out.clear();
		for (int i = 0; i < polygons.size(); i++) {
			if (isEnabled(i)) {
				points_out.push_back(polygons[i]);
			}
		}
		return points_out;
	}
	
	int addPoint(const Polygon<D, T> & P, int pointID = -1) {
		if (pointID < 0)
			pointID = polygons.size();
		polygons.push_back(P);
		polygons.back().setID(pointID);
		enabled.push_back(false);
		
		return pointID;
	}
	
	const Polygon<D, T>& operator [](int index) const {
		return polygons[index];
	}
	Polygon<D, T>& operator [](int index) {
		return polygons[index];
	}
	
	void clearHistory() {
		
	}
	void clearChanged() {
		
	}
	void markChanged() {
		
	}
	void invalidate() {
		
	}
};

#endif /* SHAPE_H_ */
