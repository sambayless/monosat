#include "DelaunayPolypartition.h"

#include "alg/Heap.h"
#include <set>
#include <algorithm>
#include "polypartition/polypartition.h"
#include <list>
template<>
DelaunayPolypartition<2,double>::DelaunayPolypartition(PolygonSet<2,double> & p):polygons(p){

}
template<>
DelaunayPolypartition<2,mpq_class>::DelaunayPolypartition(PolygonSet<2,mpq_class> & p):polygons(p){

}

template<>
void DelaunayPolypartition<2,double>::update(){
	TPPLPartition<double> pp;

	list<TPPLPoly<double>> input_polygons,result;

	for(int i =0;i<polygons.size();i++){

		if(polygons.isEnabled(i)){
			auto & polygon = polygons[i];
			input_polygons.push_back(TPPLPoly<double>());
			TPPLPoly<double> & poly = input_polygons.back();
			poly.Init(polygon.size());
			for(auto & p:polygon){
				poly[i].x = p.x;
				poly[i].y = p.y;
			}
		}
	}

	pp.Triangulate_MONO(&input_polygons,&result);

}

template<>
void DelaunayPolypartition<2,mpq_class>::update(){
	TPPLPartition<mpq_class> pp;

	list<TPPLPoly<mpq_class>> input_polygons,result;

	for(int i =0;i<polygons.size();i++){

		if(polygons.isEnabled(i)){
			auto & polygon = polygons[i];
			input_polygons.push_back(TPPLPoly<mpq_class>());
			TPPLPoly<mpq_class> & poly = input_polygons.back();
			poly.Init(polygon.size());
			for(auto & p:polygon){
				poly[i].x = p.x;
				poly[i].y = p.y;
			}
		}
	}

	pp.Triangulate_MONO(&input_polygons,&result);



}


