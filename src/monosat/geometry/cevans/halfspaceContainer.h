#ifndef HALFSPACECONTAINER_H
#define HALFSPACECONTAINER_H

#include <list>
#include <iostream>
using namespace std;

namespace cevans {
/*!
 \brief Half space and a container of integer indexes to 
 points.
 */
template<typename HS, typename PT>
class halfspaceContainer {
public:
	
	/** The half space. */
	HS halfspace;

	/** Indexes into a vector of T. */
	list<uint> index;

	/** Global points. */
	vector<PT> const & pts;

	/** Needs a half space and a reference to the points to 
	 be operated on. */
	halfspaceContainer(HS const & halfspace_, vector<PT> const & pts_) :
			halfspace(halfspace_), pts(pts_) {
	}
	
	/** Copies the targets indexes which the half-space can 
	 see, or is on the boundary. ie populates index. */
	void isInsideOrOnBoundary(list<uint> const & target);

	/** Moves the target indexes which the half-space can
	 see into index.  Includes the boundary. */
	void subtractfrom(list<uint> & target);
	
};

//---------------------------------------------------------
//  Implementation

template<typename HS, typename PT>
void halfspaceContainer<HS, PT>::subtractfrom(list<uint> & target) {
	list<uint>::iterator i = target.begin();
	for (; i != target.end(); ++i) {
		if (halfspace.isInsideOrOnBoundary(pts[*i])) {
			index.push_back(*i);
			i = target.erase(i);
		}
	}
}

template<typename HS, typename PT>
void halfspaceContainer<HS, PT>::isInsideOrOnBoundary(list<uint> const & target) {
	list<uint>::const_iterator i = target.begin();
	list<uint>::const_iterator iend = target.end();
	for (; i != iend; ++i) {
		if (halfspace.isInsideOrOnBoundary(pts[*i]))
			index.push_back(*i);
	}
}

}
;
#endif

