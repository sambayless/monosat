#ifndef INDEXTABLE_H
#define INDEXTABLE_H

#include <cassert>
#include <vector>

using namespace std;


/*!
\brief Sort a vector without altering the original vector.
An array of indexes are sorted.
\par Example
\verbatim
  int ai[10]={15,12,13,14,18,11,10,17,16,19};
  uintc aimax=10;
  vector<int> ai2(ai,ai+aimax);  // data
  vector<uint> aidxtbl2; // index array
  sortIndexTable(aidxtbl2,greater<int>(),ai2);
\endverbatim
*/
template< class Compare, class T >
void sortIndexTable
(
  vector<uint> & out,
  Compare c,
  vector<T> const & v
);


/*!
\brief Indexed containers.

I modified the code from an online C++ coding review.
Herb Sutter http://www.gotw.ca/gotw/073.htm

Changes: added a general Compare template parameter.
 Modified the variables order to my workspaces convention
 which has the output variables at the start of the function
 signature. Renamed the class variables and added an explicit
 vector interface which although in the example adds one more
 line of code I believe makes the function interface more clear.
 Chelton Evans 2007-02-06.

\par Example
\verbatim
  uintc aimax=10;
  int ai[10]={15,12,13,14,18,11,10,17,16,19};
  vector<int> aidxtbl(aimax);
  sortIndexTable(aidxtbl.begin(),less<int>(),ai,ai+aimax);
\endverbatim
*/
template< class IterOut, class Compare, class IterIn >
void sortIndexTable
( 
  IterOut out,
  Compare c,
  IterIn first, 
  IterIn last 
);

/** Randomly shuffles vector v1 with the index preserved.
    ie the original element at the i'th index
    can be found at v1[index[i]]. */
template< typename T >
void vectorshuffle
(
  vector<T> & v1,
  vector<uint> & index
);

//---------------------------------------------------------
// Implementation


template< class T, class U, class Compare >
struct comparepairfirst
{
  bool operator()
  (
    pair<T,U> const & a,
    pair<T,U> const & b 
  ) const
    { return Compare()(*a.first,*b.first); }
};

template< class IterOut, class Compare, class IterIn >
void sortIndexTable
( 
  IterOut out,
  Compare c,
  IterIn first, 
  IterIn last 
)
{
  vector<pair<IterIn,int> > s(last-first);
  uint imax(s.size());
  for (uint i=0; i<imax; ++i)
    s[i]=make_pair(first+i,i);
  sort
  (
    s.begin(),
    s.end(),
    comparepairfirst<IterIn,int,Compare>()
  );
  for (uint i=0; i<imax; ++i,++out)
    *out=s[i].second;
}


template< class Compare, class T >
void sortIndexTable
(
  vector<uint> & out,
  Compare c,
  vector<T> const & v
)
{
  out.resize(v.size());
  sortIndexTable(out.begin(),c,v.begin(),v.end());
}


template< typename T >
void vectorshuffle
(
  vector<T> & v1,
  vector<uint> & index
)
{
  uint n(v1.size());
  index.resize(n);
  for (uint i=0; i<n; ++i)
    index[i]=i;
  random_shuffle(index.begin(),index.end());

  vector<T> v2(v1.begin(),v1.end());
  for (uint i=0; i<n; ++i)
    v1[i] = v2[ index[i] ];
}



 
#endif
