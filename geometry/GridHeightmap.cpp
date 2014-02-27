/*
 * GridHeightmap.cpp
 *
 *  Created on: 2014-02-26
 *      Author: sam
 */

#include "GridHeightmap.h"

template<>
void GridHeightmap<3,double>::update(){
	for(int x = 0;x<width;x++){
		for (int y = 0;y<depth;y++){
			Point3D p = getPoint(x,y);

		}
	}
}


