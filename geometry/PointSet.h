
#ifndef POINTSET_H_
#define POINTSET_H_
#include "mtl/Vec.h"
#include "GeometryTypes.h"
using namespace Minisat;
/**
 * A dynamic point set
 */
template<unsigned int D,class T=double>
class PointSet{
	//Dimension of this shape
	vec<Point<D,T>> points;
	vec<bool> enabled;
	int sz;
public:

	int dimension(){
		return D;
	}
	int size()const{
		return sz;
	}
	int fullSize()const{
		return points.size();
	}
	bool isEnabled(int pointID){
		return enabled[pointID];
	}
	void setPointEnabled(int pointID, bool _enabled){
		if(isEnabled(pointID)==_enabled)
			return;
		enabled[pointID]=_enabled;
		if(enabled){
			sz++;
		}else{
			sz--;
		}
	}

	void getEnabledPoints(vec<Point<D>> & points_out){
		points_out.clear();
		for(int i = 0;i<points.size();i++){
			if(isEnabled(i)){
				points_out.push(points[i]);
			}
		}
	}

	int addPoint(const Point<D,T> & P){
		points.push(P);
		enabled.push(false);
		return points.size()-1;
	}
};



#endif /* SHAPE_H_ */
