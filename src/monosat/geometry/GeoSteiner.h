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

//This is a wrapper for computing exact minimal steiner trees by calling GeoSteiner.
//This assumes that geosteiner is installed either system wide or in the same directory as modsat.
//If not, then the functions will abort.

/*
 * Steiner.h
 *
 *  Created on: Jun 17, 2014
 *      Author: sam
 */

#ifndef GEOSTEINER_H_
#define GEOSTEINER_H_

#include "dgl/graph/DynamicGraph.h"
#include "dgl/SteinerTree.h"
#include "GeometryTypes.h"
#include "PointSet.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>
#include <string>
#include <sstream>

using namespace dgl;

class GeoSteiner {
public:
	
	PointSet<2> & points;
	int last_modification;
	double min_weight = 0;

	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;
	const char * euclidian_steiner_tree = "efst";
	const char * rectilinear_steiner_tree = "efst";
	const char * exact_steiner_tree = "bb";

	int INF;

	const int reportPolarity;
	std::vector<bool> in_tree;
	std::vector<int> tree_edges;

	char * tmp_file = nullptr;

public:
	
	int stats_full_updates;
	int stats_fast_updates;
	int stats_fast_failed_updates;
	int stats_skip_deletes;
	int stats_skipped_updates;
	int stats_num_skipable_deletions;
	double mod_percentage;

	double stats_full_update_time;
	double stats_fast_update_time;

	GeoSteiner(PointSet<2> & points, int _reportPolarity = 0) :
			points(points), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(0), last_history_clear(
					0), INF(0), reportPolarity(_reportPolarity) {
		
		mod_percentage = 0.2;
		stats_full_updates = 0;
		stats_fast_updates = 0;
		stats_skip_deletes = 0;
		stats_skipped_updates = 0;
		stats_full_update_time = 0;
		stats_fast_update_time = 0;
		stats_num_skipable_deletions = 0;
		stats_fast_failed_updates = 0;
		min_weight = -1;
		
	}
	~GeoSteiner() {
		
	}
	
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;
		
		if (!tmp_file) {
			tmp_file = tmpnam("points");
		}
		FILE * tmp_out = fopen(tmp_file, "w");
		INF = std::numeric_limits<int>::max();
		
		in_tree.clear();
		in_tree.resize(points.size(), false);
		
		min_weight = 0;
		tree_edges.clear();
		
		if (points.nEnabled() > 1) {
			
			for (int i = 0; i < points.size(); i++) {
				if (points.pointEnabled(i)) {
					fprintf(tmp_out, "%d %d\n", points[i].x, points[i].y);
				}
			}
		}
		fflush(tmp_out);
		std::stringstream ss;
		ss << euclidian_steiner_tree << " " << tmp_file << " | bb | grep @2";
		FILE * result = popen(ss.str().c_str(), "r");
		//Z RootZ %Gap RootLPs RootCPU RedMST
		//% @2 193579 193526.000000 0.02738 104 4.19 11.6360
		
		//note that this is only 'exact' up to double precision...
		min_weight = 0;
		if (result && fscanf(result, "%% @2 %lf\n", &min_weight) == 1) {
			
		} else {
			fclose(tmp_out);
			fclose(result);
			fprintf(stderr, "Call to GeoSteiner failed, Aborting!");
			exit(1);
		}
		fclose(tmp_out);
		fclose(result);
		/*char * line=nullptr;//hopefully big enough
		 unsigned long size=0;
		 if( getline(&line, &size, result)){
		 int size = scanf()
		 }
		 if(line)
		 free(line);
		 */
		fclose(tmp_out);
		
		assert(dbg_uptodate());
	}
	
	int weight() {
		
		update();
		
		return min_weight;
	}
	
	bool disconnected() {
		return false;		//euclidian steiner tree is never disconnected.
	}
	
	void getSteinerTree(std::vector<int> & edges) {
		edges = tree_edges;
	}
	
	bool dbg_uptodate() {
#ifdef DEBUG_GEOMETRY
		
#endif
		return true;
	}
	;
	
};

#endif

