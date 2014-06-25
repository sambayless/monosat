#include "PointSet.h"
#include "mtl/Sort.h"
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
		const vec<Point<2,double>> & points;
		double centerX;
		double centerY;
		clockwise_lt(const vec<Point<2,double>> & points,double centerX, double centerY):points(points),centerX(centerX),centerY(centerY){

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
		points_clockwise.push(i);
	}
	sort(points_clockwise,clockwise_lt(points,centerX,centerY));
	hasClockwise=true;
}
