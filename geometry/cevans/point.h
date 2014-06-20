#ifndef POINT_H
#define POINT_H

#include <cassert>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
using namespace std;

#include "typedefs.h"
namespace cevans{
// TODO -= issues: general template conflicting with same type.

//  Multiplication and Memory Management Coupled,
//  take care with writing expressions.
//
//  A good compiler should precompile the large include headers.

/*!
\brief 2D point.

The class provides both direct memory management though
 adding and modifying the point itself by {+,-,*,/}= operators,
 or the generation of temporaries throught {+,-,*,/} operators.

Reading and writing to strings and streams is also supported.
*/
template< typename T>
class point2
{
public:
  
  /** The x component/dimension. */
  T x;
  /** The y component/dimension. */
  T y;
  
  /** Set the default values to zero. */
  point2() : x(0), y(0) {}
  /** Construct a 2D point. */
  point2(T x0, T y0)  : x(x0), y(y0) {}

  /** Write out a point as a string. */
  operator stringc () const
  { 
    string s;
    { stringstream ss; ss << x; s+=ss.str(); }
    s += " ";
    { stringstream ss; ss << y; s+=ss.str(); }
    return s;
  }
  /** Interpret or read in the point. */
  void serializeInverse(stringc & s)
    {  stringstream ss(s); ss >> x >> y; }
  /** Read the point in that is in human friendly format with
      brackets and a comma from the string. */
  void serializeInverseBrackets(stringc & s);
  /** Write out the point. */
  ostream& print(ostream& os) const
    { return os << x << " " << y; }
  /** Read in a point from a stream. */
  istream & serializeInverse(istream & istr)
    { return istr >> x >> y; }


  // Arithmetic.
  /** Modify this point by adding to itself. */
  point2<T>& operator += (point2<T> const & p) 
    { x+=p.x; y+=p.y; return *this; }
  /** Create a tempory adding the two points. */
  point2<T> const operator + (point2<T> const & p) const
    { point2<T> z(*this); z+=p; return z; }
  /** Modify this point by subtracting from itself. */
  point2<T>& operator -= (point2<T> const & p) 
    { x-=p.x; y-=p.y; return *this; }
  /** Create a tempory subtracting the two points. */
  point2<T> const operator - (point2<T> const & p) const
    //{ point2<T> z(*this); z-=p; return z; } seg faulted?
    { return point2<T>(this->x-p.x,this->y-p.y); } // seg faulted too.
  /** Modify this point by multiplying to itself. */
  point2<T>& operator *= (point2<T> const & p) 
    { x*=p.x; y*=p.y; return *this; }
  /** Create a tempory multiplying the two points. */
  point2<T> const operator * (point2<T> const & p) const
    { point2<T> z(*this); z*=p; return z; }
  /** Modify this point by dividing to itself. */
  point2<T>& operator /= (point2<T> const & p) 
    { x/=p.x; y/=p.y; return *this; }
  /** Create a tempory dividing the two points. */
  point2<T> const operator / (point2<T> const & p) const
    { point2<T> z(*this); z/=p; return z; }


  /** Compare the two points for equality. */
  bool const operator == (point2<T> const & p) const 
    { return (x==p.x)&&(y==p.y); }

  /** Modify this point by adding component wise to itself. */
  template< typename W>
  point2<T>& operator += (W const & w) 
    { x+=w; y+=w; return *this; }
  /** Create a tempory adding component wise the two points. */
  template< typename W >
  point2<T> const operator + (W const & w) const
    { point2<T> z(*this); z.x+=w; z.y+=w; return z; }
  template< typename W>
  /** Modify this point by subtracting component wise to itself. */
  point2<T>& operator -= (W const & w) 
    { x-=w; y-=w; return *this; }
  /** Create a tempory subtracting component wise the two points. */
  template< typename W >
  point2<T> const operator - (W const & w) const
    { point2<T> z(*this); z.x-=w; z.y-=w; return z; }
  /** Modify this point by multiplying component wise to itself. */
  template< typename W>
  point2<T>& operator *= (W const & w) 
    { x*=w; y*=w; return *this; }
  /** Create a tempory multiplying component wise the two points. */
  template< typename W >
  point2<T> const operator * (W const & p) const
    { point2<T> z(*this); z.x*=p; z.y*=p; return z; }
  /** Modify this point by dividing component wise to itself. */
  template< typename W>
  point2<T>& operator /= (W const & w) 
    { x/=w; y/=w; return *this; }
  /** Create a tempory dividing component wise the two points. */
  template< typename W >
  point2<T> const operator / (W const & p) const
    { point2<T> z(*this); z.x/=p; z.y/=p; return z; }
  /** Not test applied to each element. */
  template< typename W >
  boolc operator != (W const & w) const
    { if (x!=w.x) return true; if (y!=w.y) return true; return false; }

  // Access.
  /** Access the elements with integers. */
  T const & operator[](uintc i) const
    { assert(i<2); if ((i%2)==0) return x; return y; }
  /** Access a reference to an element. */
  T& operator[](uintc i)
    { assert(i<2); if ((i%2)==0) return x; return y; }

  // Misc

  /** Rotate 90 degrees in a anti clockwise direction. */
  void rotate90()
    { T t=x; x=-y; y=t; }

  /** A distance measure. */
  T const dot() const 
    { return x*x + y*y; }

  /** The dot product of two 2D points. */
  T const dot( point2<T> const & w ) const
    { return x*w.x + y*w.y; }

  /** The pythagorean distance between two points. */
  T const distance() const
    { return sqrt(x*x + y*y); }

  /** Make this vector a unit vector. */
  void normalize()
    { assert(distance()!=0); *this *= 1.0/distance(); }

  /** Use this point as a ratio to sum two other numbers. */
  template< typename Z >
  Z const sumInRatio( Z const & A, Z const & B ) const
    { assert(x+y!=0); return (A*x+B*y)/(x+y); }

  /** Is this point on or inside the unit box. */
  boolc isinsideunitbox() const
  {
    if (x<(T)0)
      return false;
    if (y<(T)0)
      return false;
 
    if (x>(T)1)
      return false;
    if (y>(T)1)
      return false;

    return true;
  }

//  point2<T> & operator = (point2<T> const & B)
//    { x=B.x; y=B.y; return *this; }
};

//template< typename T >
//point2<T> & operator = (point2<T> & A, point2<T> const & B)
//{
//  A = B;
//  return A;
//}

/** Write the point out to the output stream. */
template< typename T>
ostream& operator << (ostream& os, point2<T> const & p) 
  { return p.print(os); }

/** Read in the point from the input stream. */
template< typename T >
istream & operator >> (istream & istr, point2<T> & p)
  { return p.serializeInverse(istr); }


/** Interpret p as a half space.  
  If a is below the line return true else return false. */
//template< typename T>
//bool const operator < (point2<T> const & a, point2<T> const & p)
//  { if (p.x==0) return true; return a.y < p.y/p.x*a.x; }

/** Interpret p as a half space.  
  If a is above the line return true else return false. */
template< typename T>
bool const operator > (point2<T> const & a, point2<T> const & p)
  { if (p.x==0) return false; return a.y > p.y/p.x*a.x; }



/*!
\brief 3D point.

Warning: A problem was exhibited in point3 + point3 . Unknown error.
 The only unusual thing was that both were references.
 (point3<T> const &) + (point3<T> const &)
*/
template< typename T>
class point3
{
public:

  /** The x component/dimension. */
  T x;
  /** The y component/dimension. */
  T y;
  /** The z component/dimension. */
  T z;

  /** Default is the zero vector. */
  point3()  
    {x = y = z = 0; }
  /** Construct a 3D point. */
  point3(T x0, T y0, T z0) 
    : x(x0), y(y0), z(z0) {}
  /** Construct a point from 2D. */
  point3(point2<T> const & p) : x(p.x), y(p.y), z(0) {}

  /** Write this point out to the string. */
  operator stringc () const
  { 
    string s;
    { stringstream ss; ss << x; s+=ss.str(); }
    s += " ";
    { stringstream ss; ss << y; s+=ss.str(); }
    s += " ";
    { stringstream ss; ss << z; s+=ss.str(); }
    return s;
  }

  /** Read the point in from the string. */
  void serializeInverse(stringc & s)
    { stringstream ss(s); ss >> x >> y >> z; }
  /** Read the point in that is in human friendly format with
      brackets and a comma from the string. */
  void serializeInverseBrackets(stringc & s);

  /** Write the point out to the stream. */
  ostream& print(ostream& os) const
    { return os << x << " " << y << " " << z; }
  /** Read the point in from the stream. */
  istream & serializeInverse(istream & istr)
    { return istr >> x >> y >> z; }

  /** Read an array into a vector of point3<T>. */ 
  static void readin
  (
    vector< point3<T> >& v, 
    T const * abeg, 
    T const * aend
  );

  // Arithmetic.
  /** Modify this point by adding component wise to itself. */
  point3<T>& operator += (point3<T> const & p)
    { x+=p.x; y+=p.y; z+=p.z; return *this; }
  /** Create a tempory adding the two points. */
  point3<T> const operator + (point3<T> const & p) const
    { point3<T> z(*this); z+=p; return z; } 
  /** Modify this point by subtracting component wise to itself. */
  point3<T>& operator -= (point3<T> const & p)
    { x-=p.x; y-=p.y; z-=p.z; return *this; }
  /** Create a tempory subtracting component wise the two points. */
  point3<T> const operator - (point3<T> const & p) const
    { point3<T> z(*this); z-=p; return z; }
  /** Modify this point by multiplying component wise to itself. */
  point3<T>& operator *= (point3<T> const & p)
    { x*=p.x; y*=p.y; z*=p.z; return *this; }
  /** Create a tempory multiplying component wise the two points. */
  point3<T> const operator * (point3<T> const & p) const
    { point3<T> t(*this); t*= p; return t; }
  /** Modify this point by dividing component wise to itself. */
  point3<T>& operator /= (point3<T> const & p)
    { x/=p.x; y/=p.y; z/=p.z; return *this; }
  /** Create a tempory dividing component wise the two points. */
  point3<T> const operator / (point3<T> const & p) const
    { point3<T> t(*this); t/=p; return t; }

  /** Compare the two points for equality. */
  bool const operator == (point3<T> const & p) const 
    { return (x==p.x)&&(y==p.y)&&(z==p.z); }

  /** Modify this point by adding component wise to itself. */
  template< typename W >
  point3<T>& operator += (W const & w)
    { x+=w; y+=w; z+=w; return *this; }
  /** Create a tempory adding component wise the two points. */
  template< typename W >
  point3<T> const operator + (W const & w) const
    { point3<T> t(*this); t.x+=w; t.y+=w; t.z+=w; return t; }
  /** Modify this point by subtracting component wise to itself. */
  template< typename W >
  point3<T>& operator -= (W const & w)
    { x-=w; y-=w; z-=w; return *this; }
  /** Create a tempory subtracting component wise the two points. */
  template< typename W >
  point3<T> const operator - (W const & w) const
    { point3<T> t(*this); t.x-=w; t.y-=w; t.z-=w; return t; }
  /** Modify this point by multiplying component wise to itself. */
  template< typename W >
  point3<T>& operator *= (W const & w)
    { x*=w; y*=w; z*=w; return *this; }
  /** Modify this point by multiplying component wise to itself. */
  template< typename W >
  point3<T> const operator * (W const & w) const
    { point3<T> t(*this); t.x*=w; t.y*=w; t.z*=w; return t; }
  /** Modify this point by dividing component wise to itself. */
  template< typename W >
  point3<T>& operator /= (W const & w)
    { x/=w; y/=w; z/=w; return *this; }
  /** Create a tempory dividing component wise the two points. */
  template< typename W >
  point3<T> const operator / (W const & w) const
    { point3<T> t(*this); t.x/=w; t.y/=w; t.z/=w; return t; }
  /** Not test applied to each element. */
  template< typename W >
  boolc operator != (W const & w) const
    { if (x!=w.x) return true; if (y!=w.y) return true; 
      if (z!=w.z) return true; return false; }

  // Access.
  /** Access the elements with integers. */
  T const & operator[](uintc i) const
    { uint k(i); k %= 3; if(k==0) return x; if(k==1) return y; return z; } 
  T& operator[](uintc i)
    { uint k(i); k %= 3; if(k==0) return x; if(k==1) return y; return z; }

  T const distance() const 
    { return sqrt(x*x + y*y + z*z); }

  T const distanceabs() const
    { return abs(x)+abs(y)+abs(z); }
  T const collapse() const
    { return x + y + z; }

  /** A distance measure. */
  T const dot() const 
    { return x*x + y*y + z*z; }
  /** The dot product between two points. */
  T const dot( point3<T> const & w) const
    { return x*w.x + y*w.y + z*w.z; }

  /** Make this a unit vector. */
  void normalize()
    { assert(distance()!=0); *this *= 1.0/distance(); }

  /** A common math operation to find a vector at right 
      angles to both u and v. */
  void crossproduct
  ( 
    point3<T> const & u,
    point3<T> const & v
  )
  {
    x=u.y*v.z-v.y*u.z;
    y=u.z*v.x-u.x*v.z;
    z=u.x*v.y-u.y*v.x;
  }

  /** Convert this point to 2D by throwing away the z component. */
  operator point2<T> const () const 
    { return point2<T>(x,y); }

  /** Is this point on or inside the unit box. */
  boolc isinsideunitbox() const
  {
    if (x<(T)0)
      return false;
    if (y<(T)0)
      return false;
    if (z<(T)0)
      return false;
 
    if (x>(T)1)
      return false;
    if (y>(T)1)
      return false;
    if (z>(T)1)
      return false;

    return true;
  }

};

template< typename T>
ostream& operator << (ostream& os, point3<T> const & p)
  { return p.print(os); }

template< typename T >
istream & operator >> (istream & istr, point3<T> & p)
  { return p.serializeInverse(istr); }



/*!
\brief 4D point.
*/
template< typename T>
class point4
{
public:

  T x;
  T y;
  T z;
  T w;

  point4()  
    {x = y = z = w = 0; }
  point4(T x0, T y0, T z0, T w0) 
    : x(x0), y(y0), z(z0), w(w0) {}
  point4(point4<T> const & p) : x(p.x), y(p.y), z(p.z), w(p.w) {}

  /** Write this point out to the string. */
  operator stringc () const
  { 
    string s;
    { stringstream ss; ss << x; s+=ss.str(); }
    s += " ";
    { stringstream ss; ss << y; s+=ss.str(); }
    s += " ";
    { stringstream ss; ss << z; s+=ss.str(); }
    s += " ";
    { stringstream ss; ss << w; s+=ss.str(); }
    return s;
  }
  void serializeInverse(stringc & s)
    {  stringstream ss(s); ss >> x >> y >> z >> w; }
  ostream & print(ostream & os) const
    { return os << x << " " << y << " " << z << " " << w; }
  istream & serializeInverse(istream & istr)
    { return istr >> x >> y >> z >> w; }

  // Arithmetic.
  point4<T>& operator += (point4<T> const & p)
    { x+=p.x; y+=p.y; z+=p.z; w+=p.w; return *this; }
  point4<T> operator + (point4<T> const & p) const
    { point4<T> z(*this); z+=p; return z; } 
  point4<T>& operator -= (point4<T> const & p)
    { x-=p.x; y-=p.y; z-=p.z; w-=p.w; return *this; }
  point4<T> const operator -= (point4<T> const & p) const
    { point4<T> z(*this); z-=p; return z; }
  point4<T>& operator *= (point4<T> const & p)
    { x*=p.x; y*=p.y; z*=p.z; w*=p.w; return *this; }
  point4<T> const operator * (point4<T> const & p) const
    { point4<T> t(*this); t*=p; return t; }
  point4<T>& operator /= (point4<T> const & p)
    { x/=p.x; y/=p.y; z/=p.z; w/=p.w; return *this; }
  point4<T> const operator / (point4<T> const & p) const
    { point4<T> t(*this); t/=p; return t; }

  bool const operator == (point4<T> const & p) const 
    { return (x==p.x)&&(y==p.y)&&(z==p.z)&&(w==p.w); }

  template< typename W >
  point4<T>& operator += (W const & q)
    { x+=q; y+=q; z+=q; w+=q; return *this; }
  template< typename W >
  point4<T> const operator + (W const & q) const
    { point4<T> t(*this); t.x+=q; t.y+=q; t.z+=q; t.w+=q; return t; }
//  template< typename W >
//  point4<T>& operator -= (W const & q)
//    { x-=q; y-=q; z-=q; w-=q; return *this; }
//  template< typename W >
//  point4<T> const operator - (W const & q) const
//    { point4<T> t(*this); t.x-=q; t.y-=q; t.z-=q; t.w-=q; return t; }
  template< typename W >
  point4<T>& operator *= (W const & q)
    { x*=q; y*=q; z*=q; w*=q; return *this; }
  template< typename W >
  point4<T> const operator * (W const & q) const
    { point4<T> t(*this); t.x*=q; t.y*=q; t.z*=q; t.w*=q; return t; }
  template< typename W >
  point4<T>& operator /= (W const & q)
    { x/=q; y/=q; z/=q; w/=q; return *this; }
  template< typename W >
  point4<T> const operator / (W const & q) const
    { point4<T> t(*this); t.x/=q; t.y/=q; t.z/=q; t.w/=q; return t; }
  /** Not test applied to each element. */
  template< typename W >
  boolc operator != (W const & q) const
    { if (x!=q.x) return true; if (y!=q.y) return true; 
      if (z!=q.z) return true; if (w!=q.w) return true; return false; }

  T const dot() const
    { return x*x + y*y + z*z + w*w; }
  T const dot( point4<T> const & q ) const
    { return x*q.x + y*q.y + z*q.z + w*q.w; }


  // Access.
  T const & operator[](uintc i) const
    { uint k(i); k %= 4; if(k==0) return x; 
      if(k==1) return y; if(k==2) return z; return w; } 
  T& operator[](uintc i)
    { uint k(i); k %= 4; if(k==0) return x; 
      if(k==1) return y; if(k==2) return z; return w; }

  // Conversion
  operator point3<T> const () const 
    { return point3<T>(x,y,z); }
  

};

template< typename T >
ostream& operator << (ostream& os, point4<T> const & p)
  { return p.print(os); }

template< typename T >
istream & operator >> (istream & istr, point4<T> & p)
  { return p.serializeInverse(istr); }


// --------------------------------------------------------  
// Implementation

template<class T>
void point3<T>::readin
(
  vector< point3<T> >& v, 
  T const * abeg, 
  T const * aend
)
{
  assert((aend-abeg)%3==0);

  uint n = (aend-abeg)/3;
  vector< point3<T> > v2(n);
  for (uint i=0; i<n; ++i)
  {
    v2[i] = point3<T>(*abeg,*(abeg+1),*(abeg+2));
    ++abeg; ++abeg; ++abeg;
  }

  v = v2;
}

template< typename T>
void point2<T>::serializeInverseBrackets(stringc & s)
{
  uintc sz(s.size());
  string s2(s);

  for (uint i=0; i<sz; ++i)
  {
    char ch = s[i];
    if (ch=='(')
      s2[i] = ' '; 
    if (ch==')')
      s2[i] = ' '; 
    if (ch==',')
      s2[i] = ' '; 
  }

  serializeInverse(s2);
}

template< typename T>
void point3<T>::serializeInverseBrackets(stringc & s)
{
  uintc sz(s.size());
  string s2(s);

  for (uint i=0; i<sz; ++i)
  {
    char ch = s[i];
    if (ch=='(')
      s2[i] = ' '; 
    if (ch==')')
      s2[i] = ' '; 
    if (ch==',')
      s2[i] = ' '; 
  }

  serializeInverse(s2);
}

};
#endif

