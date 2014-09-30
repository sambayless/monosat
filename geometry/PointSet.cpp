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

#include "PointSet.h"
#include <algorithm>
template<>
void PointSet<2,double>::buildClockwise(){
	points_clockwise.clear();
	double centerX = 0;
	double centerY = 0;
	for(auto & p:points){
		centerX+=p.x;
		centerY+=p.y;
	}
	centerX/=points.size();
	centerY/=points.size();

	//from http://stackoverflow.com/a/6989383

	struct clockwise_lt{
		const std::vector<Point<2,double>> & points;
		double centerX;
		double centerY;
		clockwise_lt(const std::vector<Point<2,double>> & points,double centerX, double centerY):points(points),centerX(centerX),centerY(centerY){

		}

	    bool operator () (int id_a, int id_b)
		{
	    	assert(id_a>=0);assert(id_a<points.size());
	    	assert(id_b>=0);assert(id_b<points.size());
	    	auto & a = points[id_a];
	    	auto & b = points[id_b];
			if (a[0] -centerX >= 0 && b[0] -centerX < 0)
				return true;
			if (a.x -centerX < 0 && b.x -centerX >= 0)
				return false;
			if (a.x -centerX == 0 && b.x -centerX == 0) {
				if (a.y - centerY >= 0 || b.y - centerY >= 0)
					return a.y > b.y;
				return b.y > a.y;
			}

			// compute the cross product of vectors (center -> a) x (center -> b)
			int det = (a.x -centerX) * (b.y - centerY) - (b.x -centerX) * (a.y - centerY);
			if (det < 0)
				return true;
			if (det > 0)
				return false;

			// points a and b are on the same line from the center
			// check which point is closer to the center
			int d1 = (a.x -centerX) * (a.x -centerX) + (a.y - centerY) * (a.y - centerY);
			int d2 = (b.x -centerX) * (b.x -centerX) + (b.y - centerY) * (b.y - centerY);
			return d1 > d2;
		}
	};
	for(int i =0;i<points.size();i++){
		points_clockwise.push_back(i);
	}
	std::sort(points_clockwise.begin(),points_clockwise.end() ,clockwise_lt(points,centerX,centerY));
	hasClockwise=true;
}
