#include "MonotoneDelaunay.h"

#include "alg/Heap.h"
#include <set>
#include <algorithm>
template<>
void MonotoneDelaynay<2,double>::update(){

}

template<>
void MonotoneDelaynay<2,mpq_class>::update(){

}

template<>
void  MonotoneDelaynay<2,mpq_class>::buildMonotonePolygons(){

}

template<>
void  MonotoneDelaynay<2,double>::buildMonotonePolygons(){
	monotonePolygons.clear();

	struct Edge{
		int polygonID;
		int fromPointID;

	};

	struct SortMonotoneY{
		PolygonSet<2,double> & polygons;
		bool operator()(Vertex a, Vertex b)const{
			Point<2,double> & p1 = polygons[a.polygonID][a.pointID];
			Point<2,double> & p2 = polygons[b.polygonID][b.pointID];
			if(p1.y < p2.y)
				return true;
			else if (p1.y>p2.y){
				return false;
			}else if (p1.x<p2.x){
				return true;
			}else{
				return false;
			}

		}
		SortMonotoneY(PolygonSet<2,double> & polygons):polygons(polygons){}
	};

	Heap<SortLexicographic<2,double>> q(SortMonotoneY(vertices));

	struct CmpY{
		bool operator()(const Point<2,double> & p1, const Point<2,double> & p2)const{
			if(p1.y < p2.y)
				return true;
			else if (p1.y>p2.y){
				return false;
			}else if (p1.x<p2.x){
				return true;
			}else{
				return false;
			}
		}
	};

	CmpY cmp;

	int total_points = 0;
	monotonePolygons.clear();
	sorted.clear();

	for(int polygonID =0;polygonID<polygons.size();polygonID++){
		if(polygons.isEnabled(polygonID)){
			auto & polygon = polygons[polygonID];
			total_points+=polygon.size();
			for(int pointID = 0;pointID<polygon.size();pointID++){

				//start from the vertex that is farthest to the top and left in the polygon.
				//the polygon interior must be down and to the right of that vertex

				 Point<2,double> & prev =polygon[pointID-1];
				 Point<2,double> & p = polygon[pointID];
				 Point<2,double> & next = polygon[pointID+1];
				 int prevIndex = monotonePolygons.size()-1;
				 int index =  monotonePolygons.size();
				 int nextIndex = monotonePolygons.size()+1;
				 monotonePolygons.push_back(MonotoneVertex(polygonID,pointID,prevIndex,nextIndex));

				 double det = (p[0]*prev[1] - prev[0]*p[1]);
				 VertexType type;
				 bool prevBelow = cmp(prev,p);
				 bool nextBelow = cmp(next,p);
				if(prevBelow && nextBelow){
					//a vertex is a start vertex if both neighbours are below it and the interior angle in the polygon is less than pi
					//the interior angle is less than pi if
					if(convex(next,prev,p)){
						//start vertex
						type=START;
					}else{
						//split vertex
						type=SPLIT;
					}
				}else if (!prevBelow && ! nextBelow){
					if(convex(next,prev,p)){
						//end vertex vertex
						type=END;
					}else{
						//merge vertex
						type=MERGE;
					}
				}else{
					//regular vertex
					type=REGULAR;
				}
				sorted.push_back({polygonID,pointID,index,type});
			}
		}
	}
	std::sort(sorted.begin(),sorted.end(),SortMonotoneY(polygons));

	std::vector<std::vector<Vertex>> helpers;
	helpers.resize(polygons.size());
	for(int i = 0;i<helpers.size();i++){
		helpers[i].resize(polygons[i].size());
	}

	std::set<Edge> edgeTree;



	for (auto i = sorted.rbegin(); i != sorted.rend(); ++i){
		Vertex v = sorted.back();
		sorted.pop_back();
		//check the type of the vertex
		Polygon<2,double> & polygon = polygons[v.polygonID];
		Point<2,double> & prev =polygon[v.pointID-1];
		Point<2,double> & p = polygon[v.pointID];
		Point<2,double> & next = polygon[v.pointID+1];
		int prevID = v.pointID>0 ? v.pointID-1:polygon.size()-1;
		int nextID =  (v.pointID+1)%polygon.size();


		switch(v.type){
			case START:
				{
					//Insert ei in T and set helper(ei) to vi.
					edgeTree.insert({v.polygonID,v.pointID});
					helpers[v.polygonID][v.pointID] = v;
				}
				break;
			case END:
				{
					//if helper(ei-1) is a merge vertex
					if(helpers[v.polygonID][prevID].type==MERGE) {
						//Insert the diagonal connecting vi to helper(ei-1) in D.
						addDiagonal(v.index,helpers[v.polygonID][prevID].index);
					}
					//Delete ei-1 from T
					edgeTree.erase({v.polygonID,prevID});
				}
				break;

			case SPLIT:
				//Search in T to find the edge e j directly left of vi.
				{
					Edge ei{v.polygonID,v.pointID};
					auto edgeIter = edgeTree.lower_bound(ei);
					assert(edgeIter != edgeTree.begin());

					Edge prevEdge = *(edgeIter--);
					//Insert the diagonal connecting vi to helper(ej) in D.

					int helperIndex =  helpers[prevEdge.polygonID][prevEdge.fromPointID].index;
					monotonePolygons[v.index].next =helperIndex;
					monotonePolygons[helperIndex].prev = v.index;

					helpers[prevEdge.polygonID][prevEdge.fromPointID] = v;
					helpers[v.polygonID][v.pointID] = v;
					//Insert ei in T and set helper(ei) to vi.
					edgeTree.insert(ei);

				}
				break;

			case MERGE:
				{
					//if helper(ei-1) is a merge vertex
					if(helpers[v.polygonID][prevID].type==MERGE){
						int helperIndex =  helpers[v.polygonID][prevID].index;
						monotonePolygons[v.index].next =helperIndex;
						monotonePolygons[helperIndex].prev = v.index;
					}
					//Delete ei-1 from T.
					edgeTree.erase({v.polygonID,prevID});
					//Search in T to find the edge e j directly left of vi.
					auto edgeIter = edgeTree.lower_bound({v.polygonID,v.polygonID});
					assert(edgeIter != edgeTree.begin());

					Edge prevEdge = *(edgeIter--);
					if(helpers[prevEdge.polygonID][prevEdge.fromPointID].type==MERGE){
						int helperIndex =  helpers[prevEdge.polygonID][prevEdge.fromPointID].index;
						monotonePolygons[v.index].next =helperIndex;
						monotonePolygons[helperIndex].prev = v.index;
					}

					//helper(e j)vi
					helpers[prevEdge.polygonID][prevEdge.fromPointID] = v;
				}
				break;

			case REGULAR:
				{
					//if the interior of P lies to the right of vi

						//if the interior of P lies to the right of vi
					if(Below(p,prev)) {
						//if helper(ei-1) is a merge vertex
						if(helpers[v.polygonID][prevID].type==MERGE) {
							//Insert the diagonal connecting vi to helper(ei-1) in D.
							int helperIndex =  helpers[v.polygonID][prevID].index;
							monotonePolygons[v.index].next =helperIndex;
							monotonePolygons[helperIndex].prev = v.index;
						}
						//Delete ei-1 from T.
						edgeTree.erase({v.polygonID,prevID});
						//Insert ei in T and set helper(ei) to vi.
						helpers[v.polygonID][v.pointID] = v;
						//Insert ei in T and set helper(ei) to vi.
						edgeTree.insert({v.polygonID,v.pointID});
					//SAM: this if/else condition may be incorrect, according to the text book... the text book has unclear indentation for this if/else block
					} else {
						//Search in T to find the edge ej directly left of vi.
						auto edgeIter = edgeTree.lower_bound({v.polygonID,v.polygonID});
						assert(edgeIter != edgeTree.begin());
						edgeIter--;
						Edge ej =*edgeIter;
						if(helpers[ej.polygonID][ej.fromPointID].type==MERGE){
							int helperIndex =  helpers[ej.polygonID][ej.fromPointID].index;
							monotonePolygons[v.index].next =helperIndex;
							monotonePolygons[helperIndex].prev = v.index;
						}

						//helper(e j)vi
						helpers[ej.polygonID][ej.fromPointID]=v;
					}
				}
				break;
		}
	}
	sorted.clear();

}

template<>
void  MonotoneDelaynay<2,double>::addDiagonal(int fromIndex, int toIndex){
		/*
		//Following the polypartition implementation, we are going to do this by creating two NEW vertices
		//which will then end up in separate partitioned polygons
		MonotoneVertex & from = monotonePolygons[fromIndex];
		MonotoneVertex & to = monotonePolygons[toIndex];
		from.next = to.index;
		to.prev = from.index;
		*/
		monotonePolygons[fromIndex].diagonals.push(toIndex);
		//monotonePolygons[toIndex].diagonals.push(fromIndex);//this is not needed
}

template<>
void  MonotoneDelaynay<2,double>::triangulate(){
	buildMonotonePolygons();
#ifndef NDEBUG
	for(auto & v:monotonePolygons)
		assert(!v.seen);
#endif
	for(auto & v:monotonePolygons){
		if(!v.seen){
			//find all the vertices in this monotone polygon.

		}
	}
}
