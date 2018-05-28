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

#ifndef MONOTONE_CONVEXHULL_H_
#define MONOTONE_CONVEXHULL_H_
#include "ConvexHull.h"
#include "PointSet.h"
#include <gmpxx.h>
#include "bounds/BoundingBox.h"
#include <algorithm>
#include <iostream>
#include "ConvexPolygon.h"

template<unsigned int D, class T>
class MonotoneConvexHull: public ConvexHull<D, T> {

	PointSet<D, T> & pointSet;
	NConvexPolygon<D, T> hull;
	int64_t last_modification = -1;

	int last_addition = 0;
	int last_deletion = 0;
	int history_qhead = 0;
	int last_history_clear = 0;
public:
	int64_t stats_skipped_updates = 0;
	int64_t stats_updates = 0;
	MonotoneConvexHull(PointSet<D, T> & p) :
			pointSet(p) {
		hull.setBoundingVolume(new BoundingBox<D, T, Polygon<D, T>>(hull));
	}

	void update() {

		if (pointSet.getModifications() <= last_modification) {
			//stats_skipped_updates++;
			assert(dbg_uptodate());
			return;
		}
		last_modification = pointSet.getModifications();

		if (pointSet.historyclears != last_history_clear) {
			last_history_clear = pointSet.historyclears;
			history_qhead = 0;
		} else if (last_modification > 0 && (pointSet.historyclears <= (last_history_clear + 1))) {
			//check if the newly enabled points are inside the hull, or if the newly disabled points are members of the hull. If not, then skip update.
			bool needsUpdate = false;
			for (int i = history_qhead; i < pointSet.history.size(); i++) {
				int index = pointSet.history[i].id;
				Point<D, T> & dt = pointSet[index];

				if (pointSet.history[i].addition && pointSet.pointEnabled(index)) {
					if (!hull.contains(dt, true)) {
						needsUpdate = true;
						break;
					}
				} else if (!pointSet.history[i].addition && !pointSet.pointEnabled(index)) {
					for (int i = 0; i < hull.size(); i++) {
						if (hull[i].getID() == dt.getID()) {
							//if a point is disabled, but it wasn't on the hull, then it must have been contained in the hull, so we don't need to update the hull.
							needsUpdate = true;
							break;
						}
					}
				}
			}
			history_qhead = pointSet.history.size();
			if (!needsUpdate) {
				stats_skipped_updates++;
				assert(dbg_uptodate());
				return;
			}
		}

		stats_updates++;

		if (D == 2) {
			update2d();
			return;
		}
		assert(false);
	}

	ConvexPolygon<D, T> & getHull() {
		update();
		return hull;
	}

private:

	bool dbg_uptodate() {
#ifdef DEBUG_GEOMETRY
		assert(hull.isConvex());
		for (int i = 0; i < pointSet.size(); i++) {
			Point<2, T> & p = pointSet[i];
			if (!pointSet.pointEnabled(i)) {
				
			} else {
				assert(hull.contains(p, true));
			}
		}
		/*	for (auto & p:hull){

		 assert(p == pointSet[p.getID()]);
		 assert(pointSet.pointEnabled(p.getID()));
		 }*/
#endif
		return true;
	}
	void update2d();
};

template<unsigned int D, class T>
void MonotoneConvexHull<D, T>::update2d() {

	//bool requires_update=false;

	std::vector<Point<D, T>> _points;
	pointSet.getEnabledPoints(_points);

	std::vector<Point<2, T>> & points = (std::vector<Point<2, T>> &) _points;

	std::sort(points.begin(), points.end(), SortLexicographic<2, T>());
	hull.clear();

	if (points.size() >= 3) {

		std::vector<Point<2, T>> & list = hull.getVertices();

		// Build lower hull
		for (int i = 0; i < points.size(); ++i) {
			while (list.size() >= 2 && crossDif(list[list.size() - 2], list[list.size() - 1], points[i]) <= 0)
				list.pop_back();
			list.push_back(points[i]);
		}

		// Build upper hull
		for (int i = points.size() - 2, t = list.size() + 1; i >= 0; i--) {
			while (list.size() >= t && crossDif(list[list.size() - 2], list[list.size() - 1], points[i]) <= 0)
				list.pop_back();
			list.push_back(points[i]);
		}

		assert(list.size() == 0 || (list[0] == list.back()));
		if (list.size())
			list.pop_back(); //for now, we aren't replicating the first vertex at the end of the polygon.
		/*	for(auto & p:list){
         std::cout<<p.x << "," << p.y << " ";
         }
         std::cout<<"\n";*/
		hull.reorderVertices();
	} else {
		for (auto & p : points)
			hull.addVertex(p);
	}
	hull.update();
	assert(dbg_uptodate());

}

#endif
