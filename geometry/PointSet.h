
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
	int num_enabled= 0;
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
	int nEnabled() const{
		return num_enabled;
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
			num_enabled++;
		}else{
			num_enabled--;
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
		sz++;
		return points.size()-1;
	}
    const Point<D,T>& operator [] (int index) const { return points[index]; }
    Point<D,T>&       operator [] (int index)       { return points[index]; }
	void clearHistory(){

	}
	void clearChanged(){

	}
	void markChanged(){

	}
	void invalidate(){

	}
};



#endif /* SHAPE_H_ */
