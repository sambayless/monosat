/*
 * RectangleGenerator.h
 *
 *  Created on: 2014-03-07
 *      Author: sam
 */

#ifndef CIRCLE_GENERATOR_H_
#define CIRCLE_GENERATOR_H_


#include "BitmapGenerator.h"
using namespace Minisat;

class CircleGenerator:public BitmapGenerator{

	CircleGenerator(int radius,int maxFill=0):BitmapGenerator(false){
		addParameter(radius);
		addParameter(maxFill);
	}

	virtual void generate(Bitmap & bitmap){
		int radius = getParameter(0);
		int fill = getParameter(1);

		int x = radius, y = 0;
	  	 int radiusError = 1-x;

	  	 //from http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
	  	 int x0 = radius;
	  	 int y0=radius;
	  	 bitmap.setOrigin(radius,radius);

		  while(x >= y)
		  {
			bitmap[x + x0][y + y0]=fill;
			bitmap[y + x0][ x + y0];
			bitmap[-x + x0][ y + y0];
			bitmap[-y + x0][ x + y0];
			bitmap[-x + x0][ -y + y0];
			bitmap[-y + x0][ -x + y0];
			bitmap[x + x0][ -y + y0];
			bitmap[y + x0][ -x + y0];
			y++;
			if (radiusError<0)
			{
			  radiusError += 2 * y + 1;
			} else {
			  x--;
			  radiusError+= 2 * (y - x + 1);
			}
		  }
	}


};


#endif /* RECTANGLEGENERATOR_H_ */
