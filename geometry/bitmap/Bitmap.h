/*
 * Bitmap.h
 *
 *  Created on: 2014-03-07
 *      Author: sam
 */


#ifndef BITMAP_H_
#define BITMAP_H_

#include "../GeometryTypes.h"
#include "mtl/Vec.h"
using namespace Minisat;

class Bitmap{
	int m_width;
	int m_height;
	int origin_x;
	int origin_y;
	vec<vec<int> > bitmap;



public:

	Bitmap():m_width(0),m_height(0),origin_x(0),origin_y(0){

	}

	vec<vec<int> > & getBits(){
		return bitmap;
	}

	void clear(){
		for(int i =0;i<m_width;i++){
			for(int j = 0;j<m_height;j++){
				bitmap[i][j]=-1;
			}
		}
	}

	void operator &= (Bitmap & other){
		for(int i =0;i<m_width;i++){
			for(int j = 0;j<m_height;j++){
				if(bitmap[i][j] != other[i][j]){
					bitmap[i][j] =-1;
				}
			}
		}
	}
	void operator |= (Bitmap & other){
		for(int i =0;i<m_width;i++){
			for(int j = 0;j<m_height;j++){
				if(bitmap[i][j] == -1){
					bitmap[i][j] = other[i][j];
				}else if (bitmap[i][j] != -1 && other[i][j] != -1 && bitmap[i][j] != other[i][j] ){
					assert(false);
				}
			}
		}
	}
	void growTo(int w, int h){
		assert(w>=0);
		assert(h>=0);

		if(m_height<h){
			h=m_height;
			for (int x = 0;x<m_width;x++){
				bitmap[x].growTo(m_height,-1);
			}
		}

		for (;m_width<w;m_width++){
			bitmap.push();
			bitmap[m_width].growTo(m_height,-1);
		}
	}

	bool inRange(int x, int y){
		return x>=0 && x<width() && y>=0 && y<height();
	}
	void set(int x, int y, int value){
		assert(x>=0);
		assert(y>=0);
		growTo(x+1,y+1);
		bitmap[x][y]=value;
	}
	int value(int x, int y){
		if(inRange(x,y))
			return bitmap[x][y];
		else
			return -1;
	}
	struct Row{
		int x;
		Bitmap & owner;
		Row(Bitmap & _owner, int _x):owner(_owner),x(_x){}
		int& operator[](int y)  const {owner.growTo(x,y); return owner.bitmap[x][y];}   // Here we get the value.
	};
	Row operator[](int x) { return Row(*this, x);}
	int width(){
		return m_width;
	}
	int height(){
		return m_height;
	}

	int getOriginX(){
		return origin_x;
	}
	int getOriginY(){
		return origin_y;
	}
	void setOrigin(int x, int y){
		origin_x=x;
		origin_y=y;
	}
};




#endif /* BITMAPGENERATOR_H_ */
