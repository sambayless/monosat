
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
	bool pointEnabled(int pointID){
		return  enabled[pointID];
	}
	void disablePoint(int pointID){
		setPointEnabled(pointID, false);
	}
	void enablePoint(int pointID){
		setPointEnabled(pointID, true);
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

	int addPoint(const Point<D,T> & P, int pointID){
		points.push(P);
		enabled.push(false);
		return points.size()-1;
	}

	void clearHistory(){

	}
	void clearChanged(){

	}
	void markChanged(){

	}
};



#endif /* SHAPE_H_ */
