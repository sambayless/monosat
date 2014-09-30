/*
 * TreeGenerator.h
 *
 *  Created on: 2014-03-07
 *      Author: sam
 */

#ifndef TREEGENERATOR_H_
#define TREEGENERATOR_H_

#include "BitmapGenerator.h"
#include "mtl/Rnd.h"
using namespace Monosat;

class TreeGenerator:public BitmapGenerator{

	TreeGenerator(int maxWidth,int maxHeight):BitmapGenerator(true){
		addParameter(maxWidth);
		addParameter(maxHeight);
	}

	virtual void generate(Bitmap & bitmap){
		int width = getParameter(0);
		int height = getParameter(1);
		int fill = 0;
		double seed = getRandomSeed();
		int nLeaves = irand(seed,6);

		//build a tree out of a central trunk and some circles for leaves...

	}


};



#endif /* TREEGENERATOR_H_ */
