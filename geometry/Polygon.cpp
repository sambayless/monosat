/*
 * Polygon.cpp
 *
 *  Created on: Aug 2, 2014
 *      Author: sam
 */
#include "Polygon.h"



#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>

#include <CGAL/minkowski_sum_2.h>
#include <CGAL/basic.h>
// GMP is installed. Use the GMP rational number-type.
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

//implementing these functions using the CGAL library for now




template<>
Polygon<2,mpq_class> * Polygon<2,mpq_class>::binary_union(Polygon<2,mpq_class>  * b,NPolygon<2,mpq_class>  * store){
	typedef CGAL::Cartesian<CGAL::Gmpq>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}

	if(!store){
		store = new NPolygon<2,mpq_class>();
	}else{
		store->clear();
	}
	static Polygon_with_holes_2 sum;
	sum.clear();
	if (CGAL::join (P1, P2, sum)) {
		  assert (sum.number_of_holes() == 0);
		  auto & boundary = sum.outer_boundary();
		  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
			  Point_2 & p = *iter;
			  store->addVertex(Point<2,mpq_class>(mpq_class( p[0].mpq()),mpq_class(p[1].mpq())));
		  }
	} else{
		//the polygons are disjoint.
		assert(false);
	}
	return store;
}

template<>
Polygon<2,mpq_class> * Polygon<2,mpq_class>::binary_intersect(Polygon<2,mpq_class>  * b,NPolygon<2,mpq_class>  * store){
	typedef CGAL::Cartesian<CGAL::Gmpq>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}

	if(!store){
		store = new NPolygon<2,mpq_class>();
	}else{
		store->clear();
	}

	static std::list<Polygon_with_holes_2>  res;
	res.clear();
	std::list<Polygon_with_holes_2>::const_iterator it;
	CGAL::intersection (P1, P2, std::back_inserter(res));
	assert(res.size()<=1);
	for (const Polygon_with_holes_2 & p :res) {
		  assert (p.number_of_holes() == 0);
		  auto & boundary = p.outer_boundary();
		  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
			  Point_2 & p = *iter;
			  store->addVertex(Point<2,mpq_class>(mpq_class( p[0].mpq()),mpq_class(p[1].mpq())));
		  }
	}
	return store;
}
template<>
Polygon<2,mpq_class> * Polygon<2,mpq_class>::binary_difference(Polygon<2,mpq_class>  * b,NPolygon<2,mpq_class>  * store){
	typedef CGAL::Cartesian<CGAL::Gmpq>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}

	if(!store){
		store = new NPolygon<2,mpq_class>();
	}else{
		store->clear();
	}

	static std::list<Polygon_with_holes_2>  res;
	res.clear();
	std::list<Polygon_with_holes_2>::const_iterator it;
	CGAL::difference (P1, P2, std::back_inserter(res));
	assert(res.size()<=1);
	for (const Polygon_with_holes_2 & p :res) {
		  assert (p.number_of_holes() == 0);
		  auto & boundary = p.outer_boundary();
		  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
			  Point_2 & p = *iter;
			  store->addVertex(Point<2,mpq_class>(mpq_class( p[0].mpq()),mpq_class(p[1].mpq())));
		  }
	}
	return store;
}
template<>
Polygon<2,mpq_class> * Polygon<2,mpq_class>::binary_minkowski_sum(Polygon<2,mpq_class>  * b,NPolygon<2,mpq_class>  * store){
	typedef CGAL::Cartesian<CGAL::Gmpq>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x.get_mpq_t(), p.y.get_mpq_t()));
	}

	if(!store){
		store = new NPolygon<2,mpq_class>();
	}else{
		store->clear();
	}

	Polygon_with_holes_2 sum=minkowski_sum_2 (P1, P2);

	assert (sum.number_of_holes() == 0);
	  auto & boundary = sum.outer_boundary();
	  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
		  Point_2 & p = *iter;
		  store->addVertex(Point<2,mpq_class>(mpq_class( p[0].mpq()),mpq_class(p[1].mpq())));
	  }

	return store;
}
template<>
Polygon<2,double> * Polygon<2,double>::binary_union(Polygon<2,double>  * b,NPolygon<2,double>  * store){

	typedef CGAL::Cartesian<double>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x, p.y));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x, p.y));
	}

	if(!store){
		store = new NPolygon<2,double>();
	}else{
		store->clear();
	}
	static Polygon_with_holes_2 sum;
	sum.clear();
	if (CGAL::join (P1, P2, sum)) {
		  assert (sum.number_of_holes() == 0);
		  auto & boundary = sum.outer_boundary();
		  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
			  Point_2 & p = *iter;
			  store->addVertex(Point<2,double>(p[0],p[1]));
		  }
	} else{
		//the polygons are disjoint.
		assert(false);
	}
	return store;
}
template<>
Polygon<2,double> * Polygon<2,double>::binary_intersect(Polygon<2,double>  * b,NPolygon<2,double>  * store){
	typedef CGAL::Cartesian<double>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x, p.y));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x, p.y));
	}

	if(!store){
		store = new NPolygon<2,double>();
	}else{
		store->clear();
	}

	static std::list<Polygon_with_holes_2>  res;
	res.clear();
	std::list<Polygon_with_holes_2>::const_iterator it;
	CGAL::intersection (P1, P2, std::back_inserter(res));
	assert(res.size()<=1);
	for (const Polygon_with_holes_2 & p :res) {
		  assert (p.number_of_holes() == 0);
		  auto & boundary = p.outer_boundary();
		  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
			  Point_2 & p = *iter;
			  store->addVertex(Point<2,double>(p[0],p[1]));
		  }
	}
	return store;
}
template<>
Polygon<2,double> * Polygon<2,double>::binary_difference(Polygon<2,double>  * b,NPolygon<2,double>  * store){
	typedef CGAL::Cartesian<double>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x, p.y));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x, p.y));
	}

	if(!store){
		store = new NPolygon<2,double>();
	}else{
		store->clear();
	}

	static std::list<Polygon_with_holes_2>  res;
	res.clear();
	std::list<Polygon_with_holes_2>::const_iterator it;
	CGAL::difference (P1, P2, std::back_inserter(res));
	assert(res.size()<=1);
	for (const Polygon_with_holes_2 & p :res) {
		  assert (p.number_of_holes() == 0);
		  auto & boundary = p.outer_boundary();
		  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
			  Point_2 & p = *iter;
			  store->addVertex(Point<2,double>(p[0],p[1]));
		  }
	}
	return store;
}
template<>
Polygon<2,double> * Polygon<2,double>::binary_minkowski_sum(Polygon<2,double>  * b,NPolygon<2,double>  * store){
	typedef CGAL::Cartesian<double>                Kernel;
	typedef Kernel::Point_2                             Point_2;
	typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

	static Polygon_2 P1;
	P1.clear();
	for (auto & p:*this){
		P1.push_back (Point_2 (p.x, p.y));
	}
	static Polygon_2 P2;
	P2.clear();
	for (auto & p:*b){
		P2.push_back (Point_2 (p.x, p.y));
	}

	if(!store){
		store = new NPolygon<2,double>();
	}else{
		store->clear();
	}

	Polygon_with_holes_2 sum=minkowski_sum_2 (P1, P2);

	assert (sum.number_of_holes() == 0);
	  auto & boundary = sum.outer_boundary();
	  for (auto iter = boundary.vertices_begin();iter!= boundary.vertices_end();++iter){
		  Point_2 & p = *iter;
		  store->addVertex(Point<2,double>(p[0],p[1]));
	  }

	return store;
}





