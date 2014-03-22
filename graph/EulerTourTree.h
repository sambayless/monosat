/*
 * EulerTourTree.h
 *
 *  Created on: 2014-02-07
 *      Author: sam
 */

#ifndef EULERTOURTREE_H_
#define EULERTOURTREE_H_
#include "Treap2.h"
#include <mtl/Vec.h>
using namespace Minisat;
using namespace ASoliman;
class EulerTourTree{

	Treap<int> tree;

	vec<Treap<int>::Node*> nodes;
		void cut(int v){
			Treap<int>::Node* n = tree.split(nodes[v]);


		}

		int findRoot(int v){
			return tree.findRoot(nodes[v])->getValue();
		}

		bool connected(int from, int to){
			return findRoot(from)==findRoot(to);
		}


		int link(int u, int v, int edgeID,int value) {
			  //Move both vertices to root
			  makeRoot(node);
			  makeRoot(otherNode);

			  //Create half edges and link them to each other
			  forward_edges.growTo(edgeID+1);
			  backward_edges.growTo(edgeID+1);

			  if(forward_edges[edgeID]==NULL){
				  forward_edges[edgeID] = new EulerHalfEdge();
				  backward_edges[edgeID] = new EulerHalfEdge();
			  }

			  forward_edges[edgeID]->value=value;
			  forward_edges[edgeID]->s=node;
			  forward_edges[edgeID]->t = otherNode;
			  forward_edges[edgeID]->opposite = edgeID;

			  backward_edges[edgeID]->value=value;
			  backward_edges[edgeID]->s=otherNode;
			  backward_edges[edgeID]->t = node;
			  backward_edges[edgeID]->opposite = edgeID;

			  //var st = new EulerHalfEdge(value, this, other, null, null)
			  //var ts = new EulerHalfEdge(value, other, this, null, st)
			  //st.opposite = ts
			  forward_edges[edgeID]->node = t.insert( node->node,forward_edges[edgeID],1);
			  backward_edges[edgeID]->node = t.insert(otherNode->node,backward_edges[edgeID],1);


			  //Link tours together
			  t.concat(node,otherNode);

			  //Return half edge
			  return  edgeID;
			}

};


#endif /* EULERTOURTREE_H_ */
