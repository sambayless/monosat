/*
 * GeneratorAnalyzer.h
 *
 *  Created on: 2014-03-08
 *      Author: sam
 */

#ifndef GENERATORANALYZER_H_
#define GENERATORANALYZER_H_



#include "../GeometryTypes.h"
#include "mtl/Vec.h"
using namespace Minisat;

//Given a generator as input, constructs a model of the under and over approximations over a range of parameter values
class GeneratorAnalyzer{
	Bitmap local;
	//Assumes the generator has already had its parameter space down-sampled
	void getUnderApprox_private(int parameter, BitmapGenerator & g, Bitmap & under_out){

		if(!g.isFixed(parameter)){
			int max = g.getMaximum(parameter);
			for(int j = 0;j<=max;j++){
				g.fixParameter(parameter,j);
				if(parameter<g.nParams()){
					getUnderApprox_private(parameter+1,g,under_out);
				}else{
					local.clear();
					g.generate(local);//can combine these operations and avoid the local bitmap, if needed
					under_out&=local;
				}
				g.unfixParameter(parameter);
			}
		}else{
			if(parameter<g.nParams()-1)
				getUnderApprox_private(parameter+1,g,under_out);
			else{
				local.clear();
				g.generate(local);
				under_out&=local;
			}
		}


	}
	void getOverApprox_private(int parameter, BitmapGenerator & g, Bitmap & over_out){

		if(!g.isFixed(parameter)){
			int max = g.getMaximum(parameter);
			for(int j = 0;j<=max;j++){
				g.fixParameter(parameter,j);
				if(parameter<g.nParams()){
					getOverApprox_private(parameter+1,g,over_out);
				}else{
					local.clear();
					g.generate(local);//can combine these operations and avoid the local bitmap, if needed
					over_out|=local;
				}
				g.unfixParameter(parameter);
			}
		}else{
			if(parameter<g.nParams()-1)
				getOverApprox_private(parameter+1,g,over_out);
			else{
				local.clear();
				g.generate(local);
				over_out|=local;
			}
		}
	}
public:
	GeneratorAnalyzer(){

	}

	//Assumes the generator has already had its parameter space down-sampled
	void getUnderApprox(BitmapGenerator & g, Bitmap & under_out){
		under_out.clear();
		getUnderApprox_private(0,g,under_out);
	}


	//Assumes the generator has already had its parameter space down-sampled
	void getOverApprox(BitmapGenerator & g, Bitmap & over_out){
		over_out.clear();
		getOverApprox_private(0,g,over_out);
	}

};





#endif /* GENERATORANALYZER_H_ */
