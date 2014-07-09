
#ifndef POINTSET_H_
#define POINTSET_H_
#include <vector>
#include "GeometryTypes.h"
/**
 * A dynamic point set
 */
template<unsigned int D,class T=double>
class PointSet{
	//Dimension of this shape
	std::vector<Point<D,T>> points;
	bool hasClockwise = false;
	std::vector<int> points_clockwise;
	std::vector<bool> enabled;
	int id;
public:
	long modifications = 0;
	int num_enabled= 0;

	int additions=0;
	int deletions=0;
	int historyclears=0;
	bool is_changed=0;

	void buildClockwise();
	struct PointsetChange{
		bool addition;
		int id;

		int mod;
		int prev_mod;
	};
	std::vector<PointsetChange> history;
public:
	PointSet(int id=-1):id(id){

	}
	int getID(){
		return id;
	}
	int dimension(){
		return D;
	}
	int size()const{
		return points.size();
	}
	int fullSize()const{
		return points.size();
	}
	int nEnabled() const{
		return num_enabled;
	}
	bool isEnabled(int pointID){
		assert(pointID>=0);assert(pointID<enabled.size());
		return enabled[pointID];
	}
	bool pointEnabled(int pointID){
		assert(pointID>=0);assert(pointID<enabled.size());
		return  enabled[pointID];
	}
	void disablePoint(int pointID){
		setPointEnabled(pointID, false);
	}
	void enablePoint(int pointID){
		setPointEnabled(pointID, true);
	}

	//returns indices of all points (including disabled points!) in clockwise order
	std::vector<int> & getClockwisePoints(){
		if(!hasClockwise){
			buildClockwise();
		}
		return points_clockwise;
	}

	void setPointEnabled(int pointID, bool _enabled){
		assert(pointID>=0);assert(pointID<enabled.size());
		if(isEnabled(pointID)==_enabled)
			return;
		modifications++;
		enabled[pointID]=_enabled;
		if(_enabled){
			additions=modifications;
			history.push_back({true,pointID,modifications,additions});
			num_enabled++;
		}else{
			deletions=modifications;
			history.push_back({false,pointID,modifications,deletions});
			num_enabled--;
		}
	}

	std::vector<Point<D,T>> & getEnabledPoints(std::vector<Point<D,T>> & points_out){
		points_out.clear();
		for(int i = 0;i<points.size();i++){
			if(isEnabled(i)){
				points_out.push_back(points[i]);
			}
		}
		return points_out;
	}

	int addPoint(const Point<D,T> & P){

		points.push_back(P);

		enabled.push_back(false);
		hasClockwise=false;
		modifications++;
		deletions=modifications;
		return points.size()-1;
	}

    const Point<D,T>& operator [] (int index) const { return points[index]; }
    Point<D,T>&       operator [] (int index)       { return points[index]; }

    long getModifications()const{
    	return modifications;
    }

	void clearHistory(){
		if(history.size()>10000){
			history.clear();
			historyclears++;
		}
	}

	//force a new modification
	void invalidate(){
		modifications++;
		additions=modifications;
		modifications++;
		deletions=modifications;
	}

	void markChanged(){
		is_changed=true;
	}

	bool changed(){
		return is_changed;
	}

	void clearChanged(){
		is_changed=false;
	}
};

template<>
void PointSet<2,double>::buildClockwise();

#endif /* SHAPE_H_ */
