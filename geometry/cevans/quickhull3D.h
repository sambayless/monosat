#ifndef QUICKHULL3D_H
#define QUICKHULL3D_H

#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

#include "halfspaceD3.h"


#include "halfspaceContainer.h"
#include "indextable.h"
namespace cevans{

/*!
\brief Compute the 3D convex hull of points.

O(n^2) complexity in time (worst case).

Data structure assumption: pts[0] is not used as 0 index 
 is used to mean that there is no point. 
*/
template< typename PT, typename D >
class quickhull3D
{
public:

  typedef halfspaceD3<PT,D> HS;
  typedef halfspaceContainer<HS,PT> HSC;

  /** Stack of half spaces to be processed. */
  vector< HSC* > cs;
  
  /** Global points. */
  vector<PT> const & pts;

  /** Points on the hull. */
  vector< uint > boundary;

  /** Compute the convex hull of the point pts. 
      The first element is a dummy value. */
  quickhull3D(vector<PT> const & pts_);

  /** Reset the algorithm. */
  void reset();

  /** Are there any more half-space containers to process? */
  boolc operator ! () const
    { return cs.empty()==false; } 

  /** Process a half-space container. */
  void operator ++ ();

private:

  /** Search three unique extreme points on the hull. For 
      non degenerate data ie unique xmin and xmax exist.*/
  void findMinMax(uint & q0, uint & q1, uint& q2) const;

  /** Partitions the half space by finding the furthest 
      point and creating two new halfspaces from this new 
      point.  Imagine a tetrahedron with the base as the 
      original half space (pointing inwards) and the two 
      new half spaces pointing outwards. For efficiency the 
      container h is modified. */
  void partition( HSC & h );

};


/*!
\brief Input to quickhull algorithm is shuffled.

O(nlogn) expected complexity in time.

\par Example
\verbatim
vector<pt3> pts;
//fill pts, then calculate the hull.
quickhull3Drandomized qh(pts);
\endverbatim
*/

// Changing algorithm, <TODO> re-implement later.
//template< typename PT, typename D >
//class quickhull3Drandomized
//{
//public:

  /** Global points. */
//  vector<PT> const & pts;

  /** Points on the hull. */
//  vector< uint > boundary;

  /** Compute the convex hull of the point pts. */
//  quickhull3Drandomized(vector<PT> const & pts_);
//};



//---------------------------------------------------------
//  Implementation.


template< typename PT, typename D >
void quickhull3D<PT,D>::findMinMax
(
  uint & q0, 
  uint & q1,
  uint & q2
) const
{
  q0= 1;
  q1= 1;
  q2= 1;
  uintc imax=pts.size();
  assert(imax>1);
  for (uint i=1; i<imax; ++i)
  {
    if (pts[i].x<pts[q0].x)
      q0=i;
    if (pts[i].x>pts[q1].x)
      q1=i;
    if (pts[i].y<pts[q2].y)
      q2=i;
  }
  if ((q0!=q2)&&(q1!=q2))
    return;

  q2=1;
  for (uint i=1; i<imax; ++i)
  {
    if (pts[i].y>pts[q2].y)
      q2=i;
  }
  if ((q0!=q2)&&(q1!=q2))
    return;

  q2=1;
  for (uint i=1; i<imax; ++i)
  {
    if (pts[i].z>pts[q2].z)
      q2=i;
  }
  if ((q0!=q2)&&(q1!=q2))
    return;

  q2=1;
  for (uint i=1; i<imax; ++i)
  {
    if (pts[i].z<pts[q2].z)
      q2=i;
  }
  if ((q0!=q2)&&(q1!=q2))
    return;

  assert(false);
}

template< typename PT, typename D >
void quickhull3D<PT,D>::operator ++ ()
{
  assert(cs.size()>0);

  HSC* hc = cs[cs.size()-1];
  cs.pop_back();

  partition(*hc);
  delete hc;
}

template< typename PT, typename D >
void quickhull3D<PT,D>::reset()
{
  boundary.clear();

  uintc n(pts.size());
  if (n<4)
    return;

  uint x0;
  uint x1;
  uint x2;

  findMinMax(x0,x1,x2);
  boundary.push_back(x0);
  boundary.push_back(x1);
  boundary.push_back(x2);

  assert(x0<n);
  assert(x1<n);
  assert(x2<n);

  list<uint> index;
  for (uint i=1; i<n; ++i)
  {
    if (i==x0)
      continue;
    if (i==x1)
      continue;
    if (i==x2)
      continue;

    index.push_back(i);
  }

  if (cs.empty()==false)
  {
    for ( uint i=0; i<cs.size(); ++i)
      { delete cs[i]; }

    cs.clear();
  } 

  HSC * h1 = 
    new HSC( HS(&pts[0],x0,x1,x2) ,pts );
    //new HSC( HS(pts[x0],pts[x1],pts[x2]) ,pts );
  h1->isInsideOrOnBoundary(index);
  cs.push_back(h1);

  HSC * h2 = 
    new HSC( HS(pts[x2],pts[x1],pts[x0]), pts );
  h2->isInsideOrOnBoundary(index);
  cs.push_back(h2);
}


template< typename PT, typename D >
quickhull3D<PT,D>::quickhull3D(vector<PT> const & pts_)
  : pts(pts_)
{
}

template< typename PT, typename D >
void quickhull3D<PT,D>::partition
(
  HSC & h
)
{
//cout << "partition" << endl;
//cout << SHOW(h.halfspace.p0) << " "
//     << SHOW(h.halfspace.p1) << " " 
//     << SHOW(h.halfspace.p2) << endl;

  list<uint> & target(h.index);
//cout << SHOW(print(target)) << endl;
  list<uint>::iterator i=target.begin();
  list<uint>::iterator iend=target.end();
  if (i==iend)
    return;

  list<uint>::iterator imax=i;
  ++i;
  D dmax = h.halfspace.distancefromhalfspace(pts[*imax]);
//cout << SHOW(dmax) << endl;
  D d;
  for (; i!=iend; ++i)
  {
    d=h.halfspace.distancefromhalfspace(pts[*i]);
//cout << SHOW(*i) << " " << SHOW(d) << endl;
    if( dmax < d )
    {
      dmax=d;
      imax=i;
    }
  }

  uint k(*imax);
//cout << SHOW(k) << endl;
  boundary.push_back(k);
  target.erase(imax);
//cout << SHOW(print(target)) << endl;


  HSC* h1 = new
    HSC( HS(pts[k],h.halfspace.p1,h.halfspace.p2), pts );
  //h1->isInsideOrOnBoundary(target);
  h1->subtractfrom(target);
  HSC* h2 = new
    HSC( HS(pts[k],h.halfspace.p0,h.halfspace.p1), pts );
  //h2->isInsideOrOnBoundary(target);
  h2->subtractfrom(target);
  HSC* h3 = new
    HSC( HS(pts[k],h.halfspace.p2,h.halfspace.p0), pts );
  //h3->isInsideOrOnBoundary(target);
  h3->subtractfrom(target);

  cs.push_back(h1);
  cs.push_back(h2);
  cs.push_back(h3);

//cout << SHOW(print(h1->index)) << endl;
//cout << SHOW(print(h2->index)) << endl;
//cout << SHOW(print(h3->index)) << endl;
}


/*
template< typename PT, typename D >
quickhull3Drandomized<PT,D>::quickhull3Drandomized
(
  vector<PT> const & pts_
)
  : pts(pts_)
{
  vector<uint> index;
  vector<PT> pts2(pts.begin(),pts.end());
  vectorshuffle(pts2,index);
  quickhull3D<PT,D> qh(pts2);
  uintc n=qh.boundary.size();
  boundary.resize(n);
  for (uint i=0; i<n; ++i)
    boundary[i] = index[qh.boundary[i]];
}
*/

};


#endif


