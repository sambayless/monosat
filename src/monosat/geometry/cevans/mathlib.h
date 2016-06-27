#ifndef MATHLIB_H
#define MATHLIB_H

//  Rewrite this file.
//

#include <cassert>
#include <iostream>
#include <string>
#include <vector>   
using namespace std;

#include "point.h"

#include "typedefs.h"
#include "zero.h"

#define PI 3.1415926535897932384626433
namespace cevans {

/** Mape the variable to be inside [0.0,1.0]. */
template<typename T>
void unitbound(T & val) {
	if (val < (T) 0) {
		val = (T) 0;
		return;
	}
	
	if (val > (T) 1.0)
		val = (T) 1.0;
}

/** Given the integer set 0,1,...,N-1 find vecnot such that all 
 the indexes not in vec are in vecnot. Solve for vecnot. */
void integersetdiff(vector<uint> & vecnot, vector<uint> & vec, uintc N);

/** n!/((n-k)!k!) */
void mathcombination(double & result, uintc n, uintc k);

//
// arcos is computed by subtracting arcsin from Pi/2.
//

/*
 template<class T>
 class arcsin
 {
 public:

 // y = arcsin(x) by k iterations of arcsin power series.
 T& operator()(T& y, T const & x, uintc k);
 //   Provide an error bound for iteration aswell <TODO>

 };
 */

/*!
 \brief In 3D given two vectors that make a plane find
 a third vector at right angles to this plane. */
class crossproduct {
public:
	
	/** On the right hand a is the thumb, b is the first
	 finger, w points outwards of palm. */
	template<class T>
	static void eval(T& w, T const & a, T const & b) {
		w[0] = a[1] * b[2] - b[1] * a[2];
		w[1] = a[2] * b[0] - a[0] * b[2];
		w[2] = a[0] * b[1] - a[1] * b[0];
	}
	
	/** On the right hand a is the thumb, b is the first
	 finger, w points outwards of palm. */
	template<class T>
	static void evalxyz(T& w, T const & a, T const & b) {
		w.x = a.y * b.z - b.y * a.z;
		w.y = a.z * b.x - a.x * b.z;
		w.z = a.x * b.y - a.y * b.x;
	}
	
};

/*!
 \brief Interval overlap tests.
 */
class intervalintersection {
public:
	
	/** Do the two intervals overlap? */
	template<class T>
	static boolc unordered(T const a0, T const a1, T const b0, T const b1);

	/** Find the common interval of intersection if it exists. */
	template<class T>
	static boolc unordered(T & c0, T & c1, T const a0, T const a1, T const b0, T const b1);

	/** Do the two intervals overlap? */
	template<class I>
	static boolc unordered(I const * i0, I const * i1) {
		return unordered(i0[0], i0[1], i1[0], i1[1]);
	}
	
	/** Test for overlap in two dimensions. 
	 The first row of I is the x interval, the second row
	 is the y interval.
	 */
	template<class I>
	static boolc unorderedD2(I const * box0, I const * box1) {
		if (unordered(box0, box1) == false)
			return false;
		
		return unordered(box0 + 2, box1 + 2);
	}
	
};

/*!
 \brief Solve 1D and 2D linear equations with one
 move variable than the number of equations.

 Find a point that satisfies the equations.

 For example plane to plane intersection gives two
 equations in three unknowns.
 */
template<typename T>
class solverInconsistent {
public:
	
	/** Find a solution for the equation a0*x+a1*y=c. */
	static boolc d1linearequ(T & x, T & y, T const a0, T const a1, T const c);

	/** Find a solution for the equations 
	 { a00*x+a01*y=c0, a10*x+a11*y=c1 }. */
	static boolc d2linearequ(T & x, T & y, T & z, T const a00, T const a01, T const a02, T const c0, T const a10,
			T const a11, T const a12, T const c1);
	
};

// <TODO> Is this necessary? Doc???
template<class T>
class surfaceplane {
public:
	
	void operator()(T & u, T & v, T & temp, T const & w) const;
	
};

// <TODO> Find out what this is and fix or delete.
template<typename T>
class convexity {
public:
	
	// x-axis
	T a;
	T b;
	// y-axis
	T af;
	T bf;

	// Sets 2 points which define the interval.
	convexity(T a_, T b_, T af_, T bf_) :
			a(a_), b(b_), af(af_), bf(bf_) {
	}
	convexity(); // costruct in bad state.
	
	// Given containers of values, verify the equality against the end points.
	void operator ()(bool& result, vector<T> const & x, vector<T> const & xf);

	// Samples n points between a and b, testing convexity equality.
	template<typename F>  // F is the function.
	void operator ()(bool& result,     // Is the data set convex?
			F f,              // function
			uint n,   // number of sample points
			T const a_,       // x value bound lower
			T const b_        // x value bound upper
			);
	
};

// <TODO> Find out if this is necessary.
/*
 \brief Orthogonal projection

 Client supplies orthogonal vectors v[i] and inner product.
 Calling operator solves a[i] for c.
 */
template<typename T,       // Data type
		typename V,       // vector
		typename IP,      // Inner product functional object: ()(T&,V&,V&)
		uint Dim  // Dimension 
>
class orthoproj {
public:
	
	T a[Dim];
	V v[Dim];  // Orthogonal vectors
	IP p;

	void operator()(V const & c);
	
};

/** Calculate the volume of the tetrahedron. */
void tetrahedronvolume(double & vol, point3<double> const & p0, point3<double> const & p1, point3<double> const & p2,
		point3<double> const & p3);

/** Calculate the area of the convex polygon. The points
 are ordered consecutively. */
void polygonconvexarea(double & area, vector<point2<double> > const & v);

//<TODO> Do something with the trianglearea functions.
//       Either add to triangle class or creat an area class.

void trianglearea(double & f, doublec x0, doublec y0, doublec x1, doublec y1, doublec x2, doublec y2);

/** Calcluate the triangles area given its three points. */
void trianglearea(double & f, point3<double> const & p0, point3<double> const & p1, point3<double> const & p2);

/*!
 \brief Determinant calculations.
 */
class determinant {
public:
	
	/** Calculate the determinant of a two by two matrix. */
	template<typename T>
	static T const calcD2(T const a00, T const a01, T const a10, T const a11) {
		return a00 * a11 - a01 * a10;
	}
	
	/** Calculate the determinant of a three by three matrix. */
	template<typename T>
	static T const calcD3(T const a00, T const a01, T const a02, T const a10, T const a11, T const a12, T const a20,
			T const a21, T const a22) {
		return a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20);
	}
};

template<typename T>
class solver {
public:
	
	/** Solve the 2D simultaneous equation.  All elements
	 are column vectors. */
	template<typename W>
	static boolc d2linearequ(W & x, W const & a0, W const & a1, W const & c) {
		T det = a0.x * a1.y - a1.x * a0.y;
		if (det * det <= zero<T>::val)
			return false;
		
		T detInv = (T) 1.0 / det;
		
		x.x = (c.x * a1.y - c.y * a1.x) * detInv;
		x.y = (c.y * a0.x - c.x * a0.y) * detInv;
		
		return true;
	}
	
	/** Solve the 2D simultaneous equation. */
	static boolc d2linearequ(T & x, T & y, T const a00, T const a01, T const a10, T const a11, T const c0, T const c1);
	
};

// <TODO> Move the d2/d3 matsolve into solver.

//  Solve linear equations for x0,x1.
//
//  [  a00  a01 ] [ x0 ] = [ c0 ] 
//     a10  a11     x1       c1
//  Returns false if determinant is equal to zero.
//
// All the vectors are column vectors.
template<typename T>
boolc d2matsolve(T & x, T const & a0, T const & a1, T const & c);

/** When T is a 3D point solving in one plane can fail,
 hence this routine will attempt to solve in the other two. */
template<typename T>
boolc d2matsolve3D(T & x, T const & a0, T const & a1, T const & c);

//  Solve linear equations for x0,x1.
//
//  [  a00  a01 ] [ x0 ] = [ c0 ] 
//     a10  a11     x1       c1
//  Returns false if determinant is equal to zero.
boolc d3matsolve(double & x0, double & x1, double & x2, doublec a00, doublec a01, doublec a02, doublec a10, doublec a11,
		doublec a12, doublec a20, doublec a21, doublec a22, doublec c0, doublec c1, doublec c2);

boolc d3matsolve(point3<double> & x, point3<double> const & a0,  // Column vector.
		point3<double> const & a1,  // Column vector.
		point3<double> const & a2,  // Column vector.
		point3<double> const & c);

//<TODO> Integrate lineintersection into line class.

//  Solve the intersection of two lines (a,b) and (c,d)
//  a + t0(b-a) = c + t1(d-c)
template<typename T, typename D>
boolc lineintersection(D & t0, D & t1, T const & a, T const & b, T const & c, T const & d);

/** a*t^2+b*t+c=0, solve for t. */
boolc solvequadratic(double & t0, double & t1, doublec a, doublec b, doublec c);

/*
 \brief Rotate a 2D point about the origin.

 The point is rotated in an anticlockwise direction.
 */
class transrotate2D {
	point2<double> r1;
	point2<double> r2;
public:
	
	/** The angle theta is in radians.  Rotation is anti clockwise. */
	transrotate2D(doublec theta);

	/** Rotate a 2D point. */
	void eval(point2<double> & p);

	/** Shift the point before and reverse after the rotation. */
	void eval(point2<double> & p, point2<double> const & shift);
	
};

void matrixmult(point3<double> & y, doublec * m, point3<double> const & x);

void matrixmult(point3<double> & y, point3<double> const & r0, point3<double> const & r1, point3<double> const & r2,
		point3<double> const & x);

//<TODO> *  Integrate teh lineSegments into the line class.
//       *  Remove zero and use zero<T>::test instead.

/** No division operation, return true if the two line
 segments (p1,p2) and (q1,q2) intersect. */
boolc lineSegmentIntersection(point2<double> const & p1, point2<double> const & p2, point2<double> const & q1,
		point2<double> const & q2, doublec zero);

//<TODO> Remove zero.
/** Generally returns true unless the lines are parallel.
 Only one division performed. */
boolc lineIntersection(double & tp, double & tq, point2<double> const & p1, point2<double> const & p2,
		point2<double> const & q1, point2<double> const & q2, doublec zero);

/** Only calculates both tp and tq if the line segments intersect.
 tp defines the intersection point p1 + (p2-p1)*tp, similarly tq.
 Only one division performed. */
boolc lineSegmentIntersection(double & tp, double & tq, point2<double> const & p1, point2<double> const & p2,
		point2<double> const & q1, point2<double> const & q2, doublec zero);

/*!
 \brief Circle and line intersection test.
 */
class circleLine {
public:
	
	/** Find the intersection of the line with points q0 and q1
	 with the circle at the origin. Assumes that the points
	 are 2D. */
	template<typename PT, typename PD>
	static boolc intersection2D(PT & p0, PT & p1, PD const radius, PT const & q0, PT const & q1);
	
};

//----------------------------------------------------------
//  Implementation

template<typename PT, typename PD>
boolc circleLine::intersection2D(PT & p0, PT & p1, PD const radius, PT const & q0, PT const & q1) {
	assert(zero<PD>::test((q0 - q1).dot()) == false);
	
	// Hard coded type. No templated solvequadratic.
	
	PT B(q1 - q0);
	
	auto a = B.dot();
	auto b = q0.dot(B) * 2.0;
	auto c = q0.dot() - radius * radius;
	double t0;
	double t1;
	bool res;
	res = solvequadratic(t0, t1, a, b, c);
	if (res == false)
		return false;
	
	p0 = q0 + B * t0;
	p1 = q0 + B * t1;
	
	return true;
}

// All column implementation.
template<typename T>
boolc d2matsolve(T & x, T const & a0, T const & a1, T const & c) {
	auto det = a0.x * a1.y - a1.x * a0.y;
	if (det == 0.0)
		return false;
	
	auto detInv = 1.0 / det;
	
	x.x = (c.x * a1.y - c.y * a1.x) * detInv;
	x.y = (c.y * a0.x - c.x * a0.y) * detInv;
	
	return true;
}

template<typename T>
boolc d2matsolve3D(T & x, T const & a0, T const & a1, T const & c) {
	auto det = a0.x * a1.y - a1.x * a0.y;
	// <TODO> make this test dependent on number.
	if (det != 0.0) {
		auto detInv = 1.0 / det;
		
		x.x = (c.x * a1.y - c.y * a1.x) * detInv;
		x.y = (c.y * a0.x - c.x * a0.y) * detInv;
		
		return true;
	}
	
	det = a0.x * a1.z - a1.x * a0.z;
	if (det != 0.0) {
		auto detInv = 1.0 / det;
		
		x.x = (c.x * a1.z - c.z * a1.x) * detInv;
		x.y = (c.z * a0.x - c.x * a0.z) * detInv;
		
		return true;
	}
	
	det = a0.y * a1.z - a1.y * a0.z;
	if (det != 0.0) {
		auto detInv = 1.0 / det;
		
		x.x = (c.y * a1.z - c.z * a1.y) * detInv;
		x.y = (c.z * a0.y - c.x * a0.z) * detInv;
		
		return true;
	}
	
	return false;
}

template<typename T, typename D>
boolc lineintersection(D & t0, D & t1, T const & a, T const & b, T const & c, T const & d) {
	T A(b - a);
	T B(c - d);
	T C(c - a);
	
	T x;
	bool res = d2matsolve(x, A, B, C);
	t0 = x.x;
	t1 = x.y;
	
	return res;
}

template<typename T, typename V, typename IP, uint Dim>
void orthoproj<T, V, IP, Dim>::operator()(V const & c) {
	T t;
	for (uint i = 0; i < Dim; ++i) {
		p(a[i], c, v[i]);
		p(t, v[i], v[i]);
		a[i] /= t;
	}
}

template<class T> template<typename F>
void convexity<T>::operator ()(bool& result, F f, uint n, T const a_, T const b_) {
	if (n == 0) {
		result = false;
		return;
	}
	
	a = a_;
	b = b_;
	
	f(af, a);
	f(bf, b);
	
//  cout << "af=" << af << ", bf=" << bf << endl;
	
	T one(1.0); //Unbelievabley  (j+1)/(n+1)*(b-a) was somewhere converting to int.
			
	/*
	 // Building containers and forwarding.
	 vector<T> x(n);
	 vector<T> xf(n);


	 for (uint j=0; j<n; ++j)
	 {
	 x[j] = a+(j+one)/(n+one)*(b-a);
	 f(xf[j],x[j]);
	 //cout << "x[" << j << "]=" << x[j] << ", xf[" << j << "]=" << xf[j] << endl;
	 }
	 operator()(result,x,xf);
	 */

	T x;
	T xf;
	T lambda;
	
	for (uint j = 0; j < n; ++j) {
		x = a + (j + one) / (n + one) * (b - a);
		f(xf, x);
//    cout << "x[" << j << "]=" << x << ", xf[" << j << "]=" << xf << endl;
		lambda = (x - b) / (a - b);
		if (lambda * af + bf * (one - lambda) < xf) {
			result = false;
			return;
		}
	}
	
}

// Call appropriate operator to construct the object.
template<class T>
convexity<T>::convexity() {
}

template<class T>
void convexity<T>::operator()(bool& result, vector<T> const & x, vector<T> const & xf) {
	result = true;
	
	T one(1.0);
	
	T lambda;
	
	for (uint i = 0; i < x.size(); ++i) {
		lambda = (x[i] - b) / (a - b);
		if (lambda * af + bf * (one - lambda) < xf[i]) {
			result = false;
			return;
		}
	}
}

template<class T>
void surfaceplane<T>::operator()(T & u, T & v, T & temp, T const & w) const {
	// Right shift.
	temp[1] = w[0];
	temp[2] = w[1];
	temp[0] = w[2];
	
	crossproduct::eval(u, w, temp);
	crossproduct::eval(v, u, w);
	
	//crossprod<T>()(u,w,temp);
	//crossprod<T>()(v,u,w);
}

/*
 template<class T>
 T& arcsin<T>::operator()(T& y, T const & x, uintc k)
 {
 T x2 = x * x;
 T a;
 a = 1.0; // Numerator pattern 1  1.3  1.3.5  1.3.5.7 ...
 T b;
 b = 2.0; // Denominator pattern 2  2.4  2.4.6  2.4.6.8 ...
 T xi = x; // Current power of x
 uint c = 2;

 // An improvement would be to calculate from the smallest
 // magnitude to the largest magnitude, minimising errors from
 // largely different magnitudes.

 y = x; // First term approximation.
 for (uint i=1; i<k; ++i)
 {
 xi *= x2;
 y += xi*a/(b*(c+1));
 c += 2;
 a *= (c-1);
 b *= c;
 }

 return y;
 }
 */

template<class T>
boolc intervalintersection::unordered(T const a0, T const a1, T const b0, T const b1) {
	// Order the points.
	T a[2];
	if (a0 < a1) {
		a[0] = a0;
		a[1] = a1;
	} else {
		a[1] = a0;
		a[0] = a1;
	}
	T b[2];
	if (b0 < b1) {
		b[0] = b0;
		b[1] = b1;
	} else {
		b[1] = b0;
		b[0] = b1;
	}
	
	if (a[1] < b[0])
		return false;
	
	if (b[1] < a[0])
		return false;
	
	return true;
}

template<class T>
boolc intervalintersection::unordered(T & c0, T & c1, T const a0, T const a1, T const b0, T const b1) {
	// Order the points.
	T a[2];
	if (a0 < a1) {
		a[0] = a0;
		a[1] = a1;
	} else {
		a[1] = a0;
		a[0] = a1;
	}
	T b[2];
	if (b0 < b1) {
		b[0] = b0;
		b[1] = b1;
	} else {
		b[1] = b0;
		b[0] = b1;
	}
	
//cout << SHOW(a[0]) << " " << SHOW(a[1]) << endl;
//cout << SHOW(b[0]) << " " << SHOW(b[1]) << endl;
	
	if (a[1] < b[0])
		return false;
	
	if (b[1] < a[0])
		return false;
	
	// Find the intersection interval.
	if (a[0] < b[0]) {
//cout << "a[0]<b[0]" << endl;
		if (b[1] < a[1]) {
//cout << "b[1]<a[1]" << endl;
			c0 = b[0];
			c1 = b[1];
		} else {
			c0 = b[0];
			c1 = a[1];
		}
	} else {
//cout << "b[0]<a[0]" << endl;
		if (a[1] < b[1]) {
			c0 = a[0];
			c1 = a[1];
		} else {
			c0 = a[0];
			c1 = b[1];
		}
	}
	
	return true;
}

template<typename T>
boolc solverInconsistent<T>::d1linearequ(T & x, T & y, T const a0, T const a1, T const c) {
	if (zero<T>::test(a0)) {
		x = 0;
		if (zero<T>::test(a1)) {
			y = 0;
			return (c == 0);
		}
		y = c / a1;
		return true;
	}
	
	if (zero<T>::test(a1)) {
		y = 0;
		x = c / a0;
		return true;
	}
	
	x = 1.0;
	y = (c - a0) / a1;
	return true;
}

template<typename T>
boolc solverInconsistent<T>::d2linearequ(T & x, T & y, T & z, T const a00, T const a01, T const a02, T const c0,
		T const a10, T const a11, T const a12, T const c1) {
	if (zero<T>::test(a00)) {
		if (zero<T>::test(a10)) {
			x = 0;
			return solver<T>::d2linearequ(y, z, a01, a02, a10, a11, c0, c1);
		}
		
		return d2linearequ(x, y, z, a10, a11, a12, c1, a00, a01, a02, c0);
	}
	
	T e[3];
	e[0] = a11 - a10 * a01 / a00;
	e[1] = a12 - a10 * a02 / a00;
	e[2] = c1 - a10 * c0 / a00;
	
	bool res = d1linearequ(y, z, e[0], e[1], e[2]);
	if (res == false)
		return false;
	
	x = (c0 - a01 * y - a02 * z) / a00;
	
	return true;
}

template<typename T>
boolc solver<T>::d2linearequ(T & x, T & y, T const a00, T const a01, T const a10, T const a11, T const c0, T const c1) {
	T det = a00 * a11 - a01 * a10;
	if (det * det <= zero<T>::val)
		return false;
	
	T detInv = (T) 1.0 / det;
	
	x = (c0 * a11 - c1 * a01) * detInv;
	y = (c1 * a00 - c0 * a10) * detInv;
	
	return true;
}

}
;

#endif

