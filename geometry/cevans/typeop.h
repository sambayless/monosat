#ifndef TYPEOP_H
#define TYPEOP_H
namespace cevans{
/*!
\brief Template type operator class.

At compile time templates type can be extracted or manipulated.
 This is critical in using template techniques.

For example removing a reference of a template argument to use 
 the type to access the classes typedef template declaration.

\par Example
\verbatim
template< typename FN, typename XI, typename T >
class explore 
{
public:

  typedef T Ttype;
...

template< typename EXP >
class pattern
{
...
  // Strip away the reference in EXP if it exists.
  cirbuffarr< typename typeop<EXP>::Tbare::Ttype > xi0;
...
};

  // Use of a reference causes compiler problems as a type.  
  pattern< explore<parab2 &,double*,double> & > pat(g);
\endverbatim

The idea for the classes comes from C++ Templates The Complete 
 Guide by D. Vandevoorde and N. Josuttis ISBN 0-201-73484-2.
 15.2.3 References and Qualifiers, pages 268-271.
*/
;
template< typename T >
class typeop 
{
public:

  typedef T Targ;
  typedef T Tbare;
  typedef T const Tconst;
  typedef T & Tref;
  typedef T & Tbareref;
  typedef T const & Tconstref;
  typedef T* Tptr;
  typedef T* Tbareptr;
  typedef T const * Tconstptr;
};

template< typename T >
class typeop< T & >
{
public:

  typedef T & Targ;
  typedef T Tbare;
  typedef T const Tconst;
  typedef T & Tref;
  typedef T const & Tconstref;
  typedef T* Tptr;
  typedef T* Tbareptr;
  typedef T const * Tconstptr;
};

template< typename T >
class typeop< T const > 
{
public:

  typedef T const Targ;
  typedef T Tbare;
  typedef T const Tconst;
  typedef T & Tref;
  typedef T & Tbareref;
  typedef T const & Tconstref;
  typedef T* Tptr;
  typedef T const * Tbareptr;
  typedef T const * Tconstptr;

};

template< typename T >
class typeop< T * >
{
public:

  typedef T * Targ;
  typedef T Tbare;
  typedef T const Tconst;
  typedef T & Tref;
  typedef T & Tbareref;
  typedef T const & Tconstref;
  typedef T* Tptr;
  typedef T* Tbareptr;
  typedef T const * Tconstptr;
};

template< typename T >
class typeop< T const * >
{
public:

  typedef T const * Targ;
  typedef T Tbare;
  typedef T const Tconst;
  typedef T & Tref;
  typedef T & Tbareref;
  typedef T const & Tconstref;
  typedef T* Tptr;
  typedef T* Tbareptr;
  typedef T const * Tconstptr;
};

// <TODO> Expand to include void and more pointer types.
//        How can a functions return type be extracted?


};


#endif

