
#ifndef POLYGONSET_H_
#define POLYGONSET_H_
#include <vector>
#include "GeometryTypes.h"
/**
 * A dynamic point set
 */
template<unsigned int D,class T=double>
class PolygonSet{
	//Dimension of this shape
	std::vector<Polygon<D,T>> polygons;
	std::vector<bool> enabled;

	int num_enabled= 0;


public:

	int dimension(){
		return D;
	}
	int size()const{
		return polygons.size();
	}

	int nEnabled() const{
		return num_enabled;
	}
	bool isEnabled(int polyID){
		return enabled[polyID];
	}


	void setPolygonEnabled(int polyID, bool _enabled){
		if(isEnabled(polyID)==_enabled)
			return;
		enabled[polyID]=_enabled;
		if(_enabled){
			num_enabled++;
		}else{
			num_enabled--;
		}
	}

	Point<D,T> & getPoint(int polyID, int pointID){
		return polygons[polyID][pointID];
	}

	std::vector<Polygon<D,T>> & getEnabledPolygons(std::vector<Polygon<D,T>> & points_out){
		points_out.clear();
		for(int i = 0;i<polygons.size();i++){
			if(isEnabled(i)){
				points_out.push_back(polygons[i]);
			}
		}
		return points_out;
	}

	int addPoint(const Polygon<D,T> & P, int pointID=-1){
		if(pointID<0)
			pointID=polygons.size();
		polygons.push_back(P);
		polygons.back().setID(pointID);
		enabled.push_back(false);

		return pointID;
	}

    const Polygon<D,T>& operator [] (int index) const { return polygons[index]; }
    Polygon<D,T>&       operator [] (int index)       { return polygons[index]; }

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
