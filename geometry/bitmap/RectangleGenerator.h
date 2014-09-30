/*
 * RectangleGenerator.h
 *
 *  Created on: 2014-03-07
 *      Author: sam
 */

#ifndef RECTANGLE_GENERATOR_H_
#define RECTANGLE_GENERATOR_H_


#include "BitmapGenerator.h"
using namespace Monosat;

class RectangleGenerator:public BitmapGenerator{

	RectangleGenerator(int maxWidth, int maxHeight,int maxFill=0):BitmapGenerator(false){
		addParameter(maxWidth);
		addParameter(maxHeight);
		addParameter(maxFill);
	}

	virtual void generate(Bitmap & bitmap){
		int width = getParameter(0);
		int height = getParameter(1);
		int fill = getParameter(2);
		for(int x =0;x<width;x++){
			for(int y = 0;y<height;y++){
				bitmap[x][y]=fill;
			}
		}
	}

};


#endif /* RECTANGLEGENERATOR_H_ */
