/*
 * EllipseGenerator.h
 *
 *  Created on: 2014-03-07
 *      Author: sam
 */

#ifndef ELLIPSE_GENERATOR_H_
#define ELLIPSE_GENERATOR_H_
#include "BitmapGenerator.h"
using namespace Minisat;

class EllipseGenerator:public BitmapGenerator{

	EllipseGenerator(int maxWidth,int maxHeight,int maxFill=0):BitmapGenerator(false){
		addParameter(maxWidth);
		addParameter(maxHeight);
		addParameter(maxFill);
	}

	virtual void generate(Bitmap & bitmap){
		int width = getParameter(0);
		int height = getParameter(1);
		int fill = getParameter(2);


	  	 //http://stackoverflow.com/a/10322607

	  	 bitmap.setOrigin(width/2,height/2);
	  	 int  originX = bitmap.getOriginX();
		 int  originY = bitmap.getOriginY();
	  	int hh = height * height;
	  	int ww = width * width;
	  	int hhww = hh * ww;
	  	int x0 = width;
	  	int dx = 0;

	  	// do the horizontal diameter
	  	for (int x = -width; x <= width; x++)
	  		bitmap[originX + x][ originY]=fill;

	  	// now do both halves at the same time, away from the diameter
	  	for (int y = 1; y <= height; y++)
	  	{
	  	    int x1 = x0 - (dx - 1);  // try slopes of dx - 1 or more
	  	    for ( ; x1 > 0; x1--)
	  	        if (x1*x1*hh + y*y*ww <= hhww)
	  	            break;
	  	    dx = x0 - x1;  // current approximation of the slope
	  	    x0 = x1;

	  	    for (int x = -x0; x <= x0; x++)
	  	    {
	  	    	bitmap[originX + x][ originY - y]=fill;
	  	    	bitmap[originX + x][ originY + y]=fill;
	  	    }
	  	}

	}


};


#endif /* ELLIPSEGENERATOR_H_ */
