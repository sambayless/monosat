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
#ifndef GRID_HEIGHTMAP_H
#define GRID_HEIGHTMAP_H

#include "Heightmap.h"
#include "monosat/mtl/Vec.h"
#include "PointSet.h"
//A simple heightmap composed of regularly spaced rectangular prisms
namespace Monosat {

template<unsigned int D, class T = double>
class GridHeightmap: public Heightmap<D, T> {
	int width;
	int depth;
	std::vector<std::vector<Point<D, T> > > grid;
	NPolygon<D, T> heightmesh;
	GridHeightmap(int width_, int depth_) :
			width(width_), depth(depth_) {
		
		for (int i = 0; i < width; i++) {
			grid.push();
			std::vector<Point<D, T>> &v = grid.back();
			for (int j = 0; j < depth; j++) {
				v.push( { i, j, -1 });
			}
		}
	}
	
	void setPoint(int x, int y, T height) {
		assert(x >= 0);
		assert(x < width);
		assert(y >= 0);
		assert(y < depth);
		grid[x][y]= {x,y,height};
	}

	Point<D,T> getPoint(int x, int y) {
		return grid[x][y];
	}

	void update();

	Polygon<D,T> & getHeightmesh() {
		return heightmesh;
	}

};
template<>
void GridHeightmap<3, double>::update();
}
;
#endif
