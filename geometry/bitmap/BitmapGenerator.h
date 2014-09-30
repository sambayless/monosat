/*
 * BitmapGenerator.h
 *
 *  Created on: 2014-03-07
 *      Author: sam
 */

#ifndef BITMAPGENERATOR_H_
#define BITMAPGENERATOR_H_

#include "../GeometryTypes.h"
#include "Bitmap.h"
using namespace Monosat;

class BitmapGenerator{

	int hasSeed;
	vec<int> params;
	vec<int> param_max;
	vec<bool> fixed;
public:

	virtual void generate(Bitmap & bitmap)=0;

	virtual BitmapGenerator(bool hasSeed){
		if(hasSeed)
			addParameter(1024,0);//set the random seed
	}

	virtual ~BitmapGenerator(){

	}

	bool hasRandomSeed()const{
		return hasSeed;
	}

	void setRandomSeed(int seed){
		if(hasRandomSeed()){
			setParameter(0,seed);
		}
	}

	int getRandomSeed()const{
		if(hasRandomSeed()){
			return getParameter(0);
		}else{
			return 0;
		}
	}

	int nParams()const{
		return params.size();
	}
	void setParameter(int parameter, unsigned int value){
		if(value>getMaximum(parameter)){
			value=getMaximum(parameter);
		}
		if(isFixed(parameter)){
			assert(params[parameter]==value);
		}
		params[parameter]=value;
	}
	unsigned int getMaximum(unsigned int parameter)const{
		return param_max[parameter];
	}

	unsigned int getParameter(int parameter)const{
		assert(parameter>=0);assert(parameter<params.size());
		return params[parameter];
	}
	bool isFixed(int parameter)const{
		return fixed[parameter];
	}

	void unfixParameter(int parameter){
		if(isFixed(parameter)){
			fixed[parameter]=false;
		}
	}

	//Fix a parameter at the given value (or at its current value, if unspecified)
	void fixParameter(int parameter, int atValue=-1){
		if(!isFixed(parameter)){
			if(atValue>=0){
				setParameter(parameter,atValue);
			}
			fixed[parameter]=true;
		}else{
			assert(atValue<0 || getParameter(parameter)==atValue);
		}
	}
protected:
	//Add an integer parameter to the generator.
	//All parameter values are unsigned integers in the range 0..maximum (inclusive)
	int addParameter(unsigned int maximum, unsigned int defaultValue=0){
		assert(defaultValue<=maximum);
		params.push(defaultValue);
		param_max.push(maximum);
		return params.size()-1;
	}

};




#endif /* BITMAPGENERATOR_H_ */
