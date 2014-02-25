#ifndef QUICKHULL2D_H
#define QUICKHULL2D_H

#include <cassert>
#include <vector>
#include <algorithm>
using namespace std;

#include "halfspaceD2.h"
#include "halfspaceContainer.h"
#include "indextable.h"
#include "typedefs.h"
#include "typeop.h"

/*!
\brief Compute the 2D convex hull of points.

O(n^2) complexity in time (worst case).
*/
template< typename PT, typename D >
class quickhull2D
{
  /** Search for the minimum and maximum point on the 
      x-axis. */
  void findMinMax(uint & xmin, uint & xmax) const;

  typedef halfspaceD2<PT,D> HS;
  typedef halfspaceContainer<HS,PT> HSC;

  /** Stack of half spaces to be processed. */
  vector< HSC* > cs;

  /** Loop to build the boundary with calls to partition. */
  void process();

  /** Partitions the half space by finding the furthest 
      point and creating two new halfspaces from this new 
      point.  Imagine a triangle with the base as the 
      original half space (pointing inwards) and the two 
      new half spaces pointing outwards. For efficiency 
      the container h is modified. */
  void partition( HSC & h );
public:
  
  /** Global points. */
  vector<PT> const & pts;

  /** Points on the hull. */
  vector< uint > boundary;

  /** Compute the convex hull of the point pts. */
  quickhull2D(vector<PT> const & _pts);
};




/*!
\brief Input to quickhull algorithm is shuffled.

O(nlogn) expected complexity in time.

\par Example
\verbatim
vector<pt2> pts;
//fill pts, then calculate the hull.
quickhull2Drandomized qh(pts);
\endverbatim
*/
template< typename PT, typename D >
class quickhull2Drandomized
{
public:

  /** Global points. */
  vector<PT> const & pts;

  /** Points on the hull. */
  vector< uint > boundary;

  /** Compute the convex hull of the point pts. */
  quickhull2Drandomized(vector<PT> const & pts_);
};


//---------------------------------------------------------
//  Implementation.

template< typename PT, typename D >
void quickhull2D<PT,D>::findMinMax
(
  uint & xmin, 
  uint & xmax
) const
{
  xmin = 0;
  xmax = 0;
  uintc imax=pts.size();
  assert(imax>0);
  for (uint i=1; i<imax; ++i)
  {
    if (pts[i].x<pts[xmin].x)
      xmin=i;
    if (pts[i].x>pts[xmax].x)
      xmax=i;
  }
}
    
template< typename PT, typename D >
quickhull2D<PT,D>::quickhull2D(vector<PT> const & pts_)
  : pts(pts_)
{
  uintc n(pts.size());
  if (n<3)
    return;

  uint x0;
  uint x1;
  findMinMax(x0,x1);
  boundary.push_back(x0);
  boundary.push_back(x1);

  assert(x0<n);
  assert(x1<n);

  list<uint> index;
  for (uint i=0; i<n; ++i)
  {
    if (i==x0)
      continue;
    if (i==x1)
      continue;

     index.push_back(i);
  }

  HSC * h1 = 
    new HSC( HS(pts[x0],pts[x1]) ,pts );
  h1->isInsideOrOnBoundary(index);

  HSC * h2 = 
    new HSC( HS(pts[x1],pts[x0]), pts );
  h2->isInsideOrOnBoundary(index);

  cs.push_back(h1);
  cs.push_back(h2);

  process();
}

template< typename PT, typename D >
void quickhull2D<PT,D>::process()
{

  uint sz = cs.size();
  HSC* hc;
  for ( ; sz!=0; sz=cs.size() )
  {
    hc = cs[sz-1];
    cs.pop_back();

    partition(*hc);
    delete hc;
  } 

}

template< typename PT, typename D >
void quickhull2D<PT,D>::partition
(
  HSC & h
)
{
  list<uint> & target(h.index);
  list<uint>::iterator i=target.begin();
  list<uint>::iterator iend=target.end();
  if (i==iend)
    return;

  list<uint>::iterator imax=i;
  ++i;
  D dmax = h.halfspace.distancefromhalfspace(pts[*imax]);
  D d;
  for (; i!=iend; ++i)
  {
    d=h.halfspace.distancefromhalfspace(pts[*i]);
    if( dmax < d )
    {
      dmax=d;
      imax=i;
    }
  }

  uint k(*imax);
  boundary.push_back(k);
  target.erase(imax);
  HSC* h1 = new
    HSC( HS(h.halfspace.p0,pts[k]), pts );
  h1->isInsideOrOnBoundary(target);
  cs.push_back(h1);

  HSC* h2 = new
    HSC( HS(pts[k],h.halfspace.p1), pts );
  h2->isInsideOrOnBoundary(target);
  cs.push_back(h2);
}

template< typename PT, typename D >
quickhull2Drandomized<PT,D>::quickhull2Drandomized
(
  vector<PT> const & pts_
)
  : pts(pts_)
{
  vector<uint> index;
  vector<PT> pts2(pts.begin(),pts.end());
  vectorshuffle(pts2,index);
  quickhull2D<PT,D> qh(pts2);
  uintc n=qh.boundary.size();
  boundary.resize(n);
  for (uint i=0; i<n; ++i)
    boundary[i] = index[qh.boundary[i]];
}



#endif


