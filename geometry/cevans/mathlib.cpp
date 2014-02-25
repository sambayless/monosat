
#include "mathlib.h"

#include <algorithm>
#include <cmath>
#include <sstream>
using namespace std;




void trianglearea
(
  double & f,
  doublec x0,
  doublec y0,
  doublec x1,
  doublec y1,
  doublec x2,
  doublec y2
)
{
  double b = point2<double>(x0-x1,y0-y1).distance();
  double c = point2<double>(x0-x2,y0-y2).distance();
  double dot = x1*x2+y1*y2;
  double z(b*c);
  f = 0.5*sqrt(z*z-dot*dot);
}


void trianglearea
(
  double & f,
  point3<double> const & p0,
  point3<double> const & p1,
  point3<double> const & p2
)
{
  point3<double> p;
  p.crossproduct(p1-p0,p2-p0);
  f = p.distance()*0.5;
}





/*

//  Solve linear equations for x0,x1.
//
//  [  a00  a01 ] [ x0 ] = [ c0 ] 
//     a10  a11     x1       c1
boolc d2matsolve
(
  double & x0,
  double & x1,
  doublec a00,
  doublec a01,
  doublec a10,
  doublec a11,
  doublec c0,
  doublec c1
)
{
  doublec det = a01*a10-a11*a00;
  if (det==0)
    return false;
  
  x0 = (a01*c1-a11*c0)/det;
  x1 = (a10*c0-a00*c1)/det;

  return true;
}
*/

/*
// All the vectors are column vectors.
boolc d2matsolve
(  
  point2<double> & x,
  point2<double> const & a0,  
  point2<double> const & a1,  
  point2<double> const & c
)
{
  doublec det = a0.x*a1.y - a1.x*a0.y;
  if (det==0.0)
    return false;

  doublec detInv = 1.0/det;

  x.x = (c.x*a1.y - c.y*a1.x)*detInv;
  x.y = (c.y*a0.x - c.x*a0.y)*detInv;

  return true;
}
*/


boolc d3matsolve
(
  double & x0,
  double & x1,
  double & x2,
  doublec a00,
  doublec a01,
  doublec a02,
  doublec a10,
  doublec a11,
  doublec a12,
  doublec a20,
  doublec a21,
  doublec a22,
  doublec c0,
  doublec c1,
  doublec c2
)
{
  assert(false);
  return false;
}


boolc d3matsolve
(
  point3<double> & x,
  point3<double> const & a0,  // Column vector.
  point3<double> const & a1,  // Column vector.
  point3<double> const & a2,  // Column vector.
  point3<double> const & c
)
{
  return
  d3matsolve
  (
    x.x,x.y,x.z,
    a0.x,a1.x,a2.x,
    a0.y,a1.y,a2.y,
    a0.z,a1.z,a2.z,
    c.x,c.y,c.z
  );
}




boolc solvequadratic
(
  double & t0,
  double & t1,
  doublec a,
  doublec b,
  doublec c
)
{
  doublec det = b*b - 4.0 * a * c;
  if (det<0.0)
    return false;

  double y = sqrt(det);
  double a2 = 0.5/a;

  t0 = (-b-y)*a2;
  t1 = (-b+y)*a2;

  return true;
}


transrotate2D::transrotate2D(doublec theta)
{
  doublec cos_t = cos(theta);
  doublec sin_t = sin(theta);
  r1 = point2<double>(cos_t,-sin_t);
  r2 = point2<double>(sin_t,cos_t);
}

void transrotate2D::eval(point2<double> & p)
{
  point2<double> const z(p);
  p.x = r1.dot(z);
  p.y = r2.dot(z);
}

void transrotate2D::eval
(
  point2<double> & p, 
  point2<double> const & shift
)
{
  p -= shift;
  eval(p);
  p += shift;
}

void matrixmult
( 
  point3<double> & y, 
  doublec * m, 
  point3<double> const & x 
)
{
  y.x = m[0]*x.x + m[1]*x.y + m[2]*x.z;
  y.y = m[3]*x.x + m[4]*x.y + m[5]*x.z;
  y.z = m[6]*x.x + m[7]*x.y + m[8]*x.z;
}

void matrixmult
( 
  point3<double> & y,
  point3<double> const & r0,
  point3<double> const & r1,
  point3<double> const & r2,
  point3<double> const & x
)
{
  y.x = r0.dot(x);
  y.y = r1.dot(x);
  y.z = r2.dot(x);
}


boolc lineSegmentIntersection
( 
  point2<double> const & p1,
  point2<double> const & p2,
  point2<double> const & q1,
  point2<double> const & q2,
  doublec zero
)
{
  point2<double> const a(p2-p1);
  point2<double> const aInv(-a.y,a.x);

//  point2<double> const aInv(p1.y-p2.y,p2.x-p1.x);

//  point2<double> const b(q2-q1);
//  point2<double> const bInv(b.y*-1.0,b.x);
  point2<double> const bInv(q1.y-q2.y,q2.x-q1.x);

  double alpha = a.dot(bInv);

  // Zero test on alpha.
  if (alpha + zero > 0.0)
  {
    if (alpha < zero)
      return false;
  }

  point2<double> const k(p1-q1);

  double w;

  if (alpha<0.0)
  {
    alpha *= -1.0;

    w = k.dot(bInv);
    if (w<0.0)
      return false;
    if (w>alpha)
      return false;

    w = k.dot(aInv);
    if (w<0.0)
      return false;

    if (w>alpha)
      return false;

    return true;
  }

  w = k.dot(bInv)*-1.0;
  if (w<0.0)
    return false;
  if (w>alpha)
    return false;

  w = k.dot(aInv)*-1.0;
  if (w<0.0)
    return false;
  if (w>alpha)
    return false;

  return true; 
}

boolc lineIntersection
(
  double & tp,
  double & tq,
  point2<double> const & p1,
  point2<double> const & p2,
  point2<double> const & q1,
  point2<double> const & q2,
  doublec zero
)
{
  point2<double> const a(p2-p1);
  point2<double> const aInv(-a.y,a.x);
  point2<double> const b(q2-q1);
  point2<double> const bInv(-b.y,b.x);

  doublec alpha = b.dot(aInv);

  // Zero test on alpha.
  if (alpha + zero > 0.0)
  {
    if (alpha < zero)
      return false;
  }

  doublec c = 1.0/alpha;

  point2<double> const w(p1-q1);
  tq = c*w.dot(aInv);
  tp = c*w.dot(bInv);

  return true;
}




boolc lineSegmentIntersection
(
  double & tp,
  double & tq,
  point2<double> const & p1,
  point2<double> const & p2,
  point2<double> const & q1,
  point2<double> const & q2,
  doublec zero
)
{
  point2<double> const a(p2-p1);
  point2<double> const aInv(-a.y,a.x);
  point2<double> const b(q2-q1);
  point2<double> const bInv(-b.y,b.x);

  doublec alpha = b.dot(aInv);

  // Zero test on alpha.
  if (alpha + zero > 0.0)
  {
    if (alpha < zero)
      return false;
  }

  doublec c = 1.0/alpha;

  point2<double> const w(p1-q1);
  tq = c*w.dot(aInv);

  if (tq<0.0)
    return false;

  if (tq>1.0)
    return false;

  tp = c*w.dot(bInv);

  if (tp<0.0)
    return false;

  if (tp>1.0)
    return false;

  return true;
}

void tetrahedronvolume
( 
  double & vol, 
  point3<double> const & p0,
  point3<double> const & p1,
  point3<double> const & p2,
  point3<double> const & p3
)
{
  //point3<double> a(p1-p0);
  //point3<double> b(p2-p0);
  //point3<double> c(p3-p0);

/*
cout << SHOW(p0) << endl;
cout << SHOW(p1) << endl;
cout << SHOW(p2) << endl;
cout << SHOW(p3) << endl;

cout << SHOW(a) << endl;
cout << SHOW(b) << endl;
cout << SHOW(c) << endl;
*/

  point3<double> q;
  crossproduct::eval(q,p1-p0,p2-p0);
  //crossprod< point3<double> >()(q,p1-p0,p2-p0);
//cout << SHOW(q) << endl;
  vol = (p3-p0).dot(q)/6.0;
}




void polygonconvexarea
(
  double & area,
  vector< point2<double> > const & v
)
{
  uintc n = v.size();

  area = v[n-1].x*v[0].y - v[0].x*v[n-1].y;
  for (uint i=0; i<n-1; ++i)
  {
    area += v[i].x*v[i+1].y - v[i+1].x*v[i].y;
  }
  area *= 0.5;
}




void mathcombination
(
  double & result,
  uintc n,
  uintc k
)
{
// (n,k) = prod(i=0..n-k-1, (n-i)/(n-k-i) )
  result = 1.0;

  uintc ndiffk = n-k;
  for (uint i=0; i<ndiffk; ++i)
  {
    result *= (n-i);
    result /= (ndiffk-i);
  }
}



void integersetdiff
(
  vector<uint> & vecnot,
  vector<uint> & vec,
  uintc N
)
{
  vecnot.clear();
  sort(vec.begin(),vec.end());

  uint k=0;
  for (uint i=0; i<N; ++i)
  {
    if (vec[k]==i)
      ++k;
    else
      vecnot.push_back(i);
  }
}








    





