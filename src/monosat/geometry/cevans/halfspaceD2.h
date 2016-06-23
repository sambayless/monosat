#ifndef HALFSPACED2_H
#define HALFSPACED2_H

#include <cassert>
#include <sstream>
using namespace std;

#include "typedefs.h"
#include "partitionspace.h"
#include "mathlib.h"
#include "zero.h"
namespace cevans {
/*!
 <TODO> Rewrite half-space storing only one point p0 and
 the normal.

 \brief Define a 2D half space with two ordered points.
 The half space points to the left of the directed line(p0,p1).

 Need to define zero. eg
 \verbatim
 template<>
 double zero<double>::val = 1E-15;
 \endverbatim
 */
template<typename PT, typename PD>
class halfspaceD2: public partitionspace<PT> {
public:
	
	/** Start point. */
	PT p0;
	/** End point. */
	PT p1;
	/** Normal but not normalized. */
	PT normal;

	/** Construct in uninitialized state. */
	halfspaceD2() {
	}
	
	/** Construct a half space from the ordered points.
	 The normal is the line rotated 90 degrees anti 
	 clockwise. */
	halfspaceD2(PT const & p0_, PT const & p1_) {
		set(p0_, p1_);
	}
	/** Construct a half space from the ordered points.
	 Left of the directed line. */
	void set(PT const & p0_, PT const & p1_);

// Problems with template, tried explicit copy constructor.
//  halfspaceD2( halfspaceD2<PT,PD> const & h)
//    : p0(h.p0), p1(h.p1), normal(h.normal) {}
	
	/** Calculate the normal the ordered points. */
	void normalcalculate() {
		PT t(p1 - p0);
		normal.x = -t.y;
		normal.y = t.x;
	}
	
	/** Is the point inside the half space? */
	boolc isInside(PT const & x) const {
		return 0 < (x.x - p0.x) * normal.x + (x.y - p0.y) * normal.y;
	}
	
	/** Is the point on or inside the half space? */
	boolc isInsideOrOnBoundary(PT const & x) const {
		return 0 < zero<PD>::val + (x.x - p0.x) * normal.x + (x.y - p0.y) * normal.y;
	}
	//{ return 0 < zero+(x.x-p0.x)*normal.x + (x.y-p0.y)*normal.y; }
	
	/** Line [p0,p1] for t=[0,1] */
	template<typename U>
	void pointOnLine(PT & x, U const t) const {
		x = p0 + (p1 - p0) * t;
	}
	
	/** Solve (t.x,t.y) : a+(-a+b)t.x = pointOnLine(t.y) */
	boolc intersectionIn_t(PT & t, PT const & a, PT const & b) const;

	/** Find the intersection of the line through (a,b) and 
	 this line. */
	boolc intersection(PT & p, PT const & a, PT const & b) const;

	/** The client configures zero for the isOnBoundary 
	 routine. */
	boolc isOnBoundary(PT const & w) const {
		PD x = (p1.y - p0.y) * (w.x - p0.x) - (p1.x - p0.x) * (w.y - p0.y);
		if (0 < (x + zero<PD>::val)) {
			if (x < zero<PD>::val)
				return true;
		}
		
		return false;
	}
	
	/** Clip the line segment if cutting the half-space. A false
	 result rejects the line as it is outside the half-space.  
	 Line a+m*t with [t0,t1]. */
	boolc clip(PD & t0, PD & t1, PT const & a, PT const & m) const;

	/** Invert the clip. */
	boolc clipNeg(PD & t0, PD & t1, PT const & a, PT const & m) const;

	/** Find the nearest point on the line to w. */
	void minimizepointtoline(PD & t, PT const & w) const {
		PT z(p1 - p0);
		t = ((w.x - p0.x) * z.x + (w.y - p0.y) * z.y) / (z.x * z.x + z.y * z.y);
	}
	
	/** A measure of the minimum distance from the line to the
	 point w without a square root. */
	PD const distancefromhalfspace(PT const & w) const;

	/** Serialize this object by writing it out as a string. */
	operator stringc() const {
		stringstream ss;
		ss << p0 << " " << p1;
		return ss.str();
	}
	
	/** For debug print this object. */
	stringc print() const {
		stringstream ss;
		ss << SHOW(p0) << " " << SHOW(p1) << " " << SHOW(normal);
		return ss.str();
	}
	
};

//---------------------------------------------------------
//  Implementation.

template<typename PT, typename PD>
PD const halfspaceD2<PT, PD>::distancefromhalfspace(PT const & w) const {
	PD t;
	minimizepointtoline(t, w);
	PT w2;
	pointOnLine(w2, t);
	PT w3(w2.x - w.x, w2.y - w.y);
	return w3.x * w3.x + w3.y * w3.y;
}

template<typename PT, typename PD>
void halfspaceD2<PT, PD>::set(PT const & p0_, PT const & p1_) {
	p0 = p0_;
	p1 = p1_;
	normalcalculate();
}

template<typename PT, typename PD>
boolc halfspaceD2<PT, PD>::intersection(PT & p, PT const & a, PT const & b) const {
	PT const u(b - a);
	PT const v(p0 - p1);
	PT const w(p0 - a);
	// t0*u+t1*v=w
	PT t;
	bool res = d2matsolve(t, u, v, w);
	if (res == false)
		return false;
	
	pointOnLine(p, t.y);
	
	return true;
}

template<typename PT, typename PD>
boolc halfspaceD2<PT, PD>::intersectionIn_t(PT & t, PT const & a, PT const & b) const {
	PT const u(b - a);
	PT const v(p0 - p1);
	PT const w(p0 - a);
	// t0*u+t1*v=w
	
	return d2matsolve(t, u, v, w);
}

template<typename PT, typename PD>
boolc halfspaceD2<PT, PD>::clipNeg(PD & t0, PD & t1, PT const & a, PT const & m) const {
	
	int k = 0;
	if (!isInside(a + m * t0))
		k += 1;
	if (!isInside(a + m * t1))
		k += 2;
	
	switch (k) {
		case 0:
			return false;
		case 1: {
			point2<PD> t;
			if (solver<PD>::d2linearequ(t, m, p0 - p1, p0 - a) == false)
				return !isInside(a);
			t1 = t[0];
			break;
		}
		case 2: {
			point2<PD> t;
			if (solver<PD>::d2linearequ(t, m, p0 - p1, p0 - a) == false)
				return !isInside(a);
			t0 = t[0];
			break;
		}
		case 3:
			return true;
	}
	
	return true;
}

template<typename PT, typename PD>
boolc halfspaceD2<PT, PD>::clip(PD & t0, PD & t1, PT const & a, PT const & m) const {
	int k = 0;
	if (isInside(a + m * t0))
		k += 1;
	if (isInside(a + m * t1))
		k += 2;
	
	switch (k) {
		case 0:
			return false;
		case 1: {
			point2<PD> t;
			if (solver<PD>::d2linearequ(t, m, p0 - p1, p0 - a) == false)
				return isInside(a);
			t1 = t[0];
			break;
		}
		case 2: {
			point2<PD> t;
			if (solver<PD>::d2linearequ(t, m, p0 - p1, p0 - a) == false)
				return isInside(a);
			t0 = t[0];
			break;
		}
		case 3:
			return true;
	}
	
	return true;
}

}
;

#endif

