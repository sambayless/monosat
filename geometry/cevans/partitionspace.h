#ifndef PARTITIONSPACE_H
#define PARTITIONSPACE_H

#include <cassert>
using namespace std;

#include "typedefs.h"

/*!
\brief A partition of a continuous n-dimension space.

This is a base class for a partition of (in general) a 
 continuous space. The template parameter T is the 
 element type in the space.

A continuous space implies that the partition has a 
 boundary.  There are two basic questions about the 
 boundary which are really each others inverses.  
 1. Given two points on either side of the partition 
    how do we calculate a point on the boundary. 
    See the intersection routine.
 2. Given a point is it on the boundary?  
    If these questions are not needed the routines do 
    not need to be overridden.
*/
template< typename PT >
class partitionspace
{
public:

  /** Is the point inside the partition? */
  virtual boolc isInside(PT const & x) const = 0; 

  /** Find the point on the boundary of the partition 
      between two points inside and outside the 
      partition. */
  virtual boolc intersection
  (
    PT & x, 
    PT const & poutside, 
    PT const & pinside
  ) const
    { assert(false); return false; }

  /** Is the point on the partition boundary? */
  virtual boolc isOnBoundary(PT const & w) const
    { assert(false); return false; }

  /** Destructor. */
  virtual ~partitionspace() {}

  enum ptclassification { undefined, below, on, above };
  /** Is the point above, below or on the boundary? */
  ptclassification classify(PT const & w) const
  { 
    if (this->isOnBoundary(w))
      return on;

    if (this->isInside(w))
      return above;

    return below;
  }
  /** Get the classification as a string. */
  static stringc classifystring(ptclassification c)
  {
    string s;
    switch (c)
    {
      case undefined: s="undefined"; break;
      case below: s="below"; break;
      case on: s="on"; break;
      case above: s="above"; break;
    }

    return s;
  } 

  /** Classify a container of points. */
  void classify
  (
    ptclassification* vc, 
    PT const* beg,
    PT const* end
  ) const;

};

//---------------------------------------------------------
// Implementation

template< typename PT >
void partitionspace<PT>::classify
(
  ptclassification* vc, 
  PT const* beg,
  PT const* end
) const
{
  PT const* i(beg);
  for (; i!=end; ++i)
  {
    *vc = classify(*i);
    ++vc;
  }
}


#endif


