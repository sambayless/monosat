/*
 * Downsampler.h
 *
 *  Created on: 2014-03-08
 *      Author: sam
 */

#ifndef DOWNSAMPLER_H_
#define DOWNSAMPLER_H_


#include "../GeometryTypes.h"
#include "Bitmap.h"
#include "mtl/Rnd.h"
#include <algorithm>
#include <cmath>
using namespace Minisat;

//Takes a high-dimensional generator and creates a version with a down-sampled and normalized parameter space
class BinarySampler:public BitmapGenerator{

	//Zudzik's pairing function:
	unsigned int pair(unsigned int a,unsigned int b){
		return a >= b ? a * a + a + b : a + b * b;
	}

	void unpair(unsigned int z,unsigned int & a_out,unsigned int & b_out){
		unsigned int root_z = (unsigned int) sqrt(z);
		unsigned int rootSquare = root_z*root_z;
		if(z-rootSquare<root_z){
			a_out = z-rootSquare;
			b_out=root_z;
		}else{
			a_out = root_z;
			b_out = z-rootSquare-root_z;
		}
	}

	unsigned long pair(unsigned long a,unsigned long b){
		return a >= b ? a * a + a + b : a + b * b;
	}

	void unpair(unsigned long z,unsigned long & a_out,unsigned long & b_out){
		unsigned long root_z = (unsigned long) sqrt(z);
		unsigned long rootSquare = root_z*root_z;
		if(z-rootSquare<root_z){
			a_out = z-rootSquare;
			b_out=root_z;
		}else{
			a_out = root_z;
			b_out = z-rootSquare-root_z;
		}
	}
public:
	BitmapGenerator & g;
	vec<int> bit_mapping;
	void generate(Bitmap & bitmap){
		int j =0;
		int p=0;
		for(int i = 0;i<bit_mapping.size();i++){
			int bits = bit_mapping[i];
			int val = 0;
			int start = j;
			for(;j<bits;j++){
				val += getParameter(j) << (j-start);
			}
			while(g.isFixed(p))
				p++;
			assert(p<g.nParams());
			g.setParameter(p,val);
		}
		g.generate(bitmap);
	}

	//Combine the parameter space of the generator into a vector of binary options
	BinarySampler(BitmapGenerator & _g):g(_g){
		int total_bits=0;

		for(int p = 0;p<g.nParams();p++){

			if(g.isFixed(p)){
				//for now, I'm going to assume that the exact same parameters will still be fixed in g later.
			}else{
				int required_bits = 0;
				int max = g.getMaximum(p);
				while (max >>= 1)
					++required_bits;
				total_bits+=required_bits;
				bit_mapping.push(total_bits);
			}
		}
		for(int i =0;i<total_bits;i++){
			addParameter(1);
		}
	}

};



#endif /* DOWNSAMPLER_H_ */
