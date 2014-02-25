#ifndef D3HALFSPACE_H
#define D3HALFSPACE_H

#include <cassert>
using namespace std;


#include "partitionspace.h"
#include "mathlib.h"


/*!
\brief Define a 3D half space with three ordered points.

The inside of the half space is on the anti clockwise
 side of the plane with the three points.
*/
template< typename PT, typename PD >
class halfspaceD3 : public partitionspace<PT>
{
public:

  /** Start point. */
  PT p0;  
  /** Second point. */
  PT p1;
  /** End point. */
  PT p2;  
  /** Normal. */
  PT normal;  

  /** Construct in uninitialized state. */
  halfspaceD3() {}

  template< typename INDX >
  halfspaceD3
  (
    PT const * pts,
    INDX i0,
    INDX i1,
    INDX i2
  )
    : p0(pts[i0]), p1(pts[i1]), p2(pts[i2])
    { normalcalculate(); }

  /** Construct the plane with anticlockwise point 
      ordering.  The normal is outwards on the 
      anticlockwise side of the plane. */
  halfspaceD3
  ( 
    PT const & p0_,
    PT const & p1_,
    PT const & p2_
  );
  /** Construct a half space from the ordered points. */
  void set
  ( 
    PT const & p0_,
    PT const & p1_,
    PT const & p2_
  );

  /** Calculate the normal the ordered points. */
  void normalcalculate();

  /** Is the point inside the half space? */
  boolc isInside( PT const & x ) const
    { return 0 < normal.x*(x.x-p0.x)+normal.y*(x.y-p0.y)+normal.z*(x.z-p0.z); }
  /** Is the point on or inside the half space? */
  boolc isInsideOrOnBoundary( PT const & x ) const
    { return 0 < zero<PD>::val+normal.x*(x.x-p0.x)+normal.y*(x.y-p0.y)+normal.z*(x.z-p0.z); }
    //{ return 0 < zero+normal.x*(x.x-p0.x)+normal.y*(x.y-p0.y)+normal.z*(x.z-p0.z); }


/*
// <TODO> Remove crossproduct.  Check what and where this func works.
  // Calculate the area of the triangle a,b,c.
  T const trianglearea() const
  {
    T c;
    c.crossproduct(p1-p0,p2-p0);
//cout << SHOW(c) << endl;

    return c.distance()*.5;
  }
*/

  /** A measure of the distance from the plane to the
      point w without a square root. ie
      N^2*d^2 where d is the vanilla distance function. */
  PD const distancefromhalfspace(PT const & w) const
    { 
      PT p;
      p.x=normal.x*normal.x*(p0.x-w.x);
      p.y=normal.y*normal.y*(p0.y-w.y);
      p.z=normal.z*normal.z*(p0.z-w.z);
      return p.x*p.x+p.y*p.y+p.z*p.z;
/*
      return normal.x*normal.x*(p0.x-w.x)+
             normal.y*normal.y*(p0.y-w.y)+
             normal.z*normal.z*(p0.z-w.z);
*/
    }


};


// --------------------------------------------------------- 
// Implementation


/*
  template< typename INDX >
  halfspaceD3
  (
    PT const * pts,
    INDX i0,
    INDX i1,
    INDX i2
  );
*/

template< typename PT, typename PD >
void halfspaceD3<PT,PD>::normalcalculate()
{
  PT u(p2 - p1);
  PT v(p0 - p1);

  normal.x=u.y*v.z-v.y*u.z;
  normal.y=u.z*v.x-u.x*v.z;
  normal.z=u.x*v.y-u.y*v.x;
}

template< typename PT, typename PD >
halfspaceD3<PT,PD>::halfspaceD3
(
  PT const & p0_,
  PT const & p1_,
  PT const & p2_
)
  : p0(p0_), p1(p1_), p2(p2_)
{
  normalcalculate();
}

template< typename PT, typename PD >
void halfspaceD3<PT,PD>::set
(
  PT const & p0_,
  PT const & p1_,
  PT const & p2_
)
{
  p0=p0_;
  p1=p1_;
  p2=p2_;

  normalcalculate();
}


#endif


