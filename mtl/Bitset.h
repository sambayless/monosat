/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/
//Barebones bit-vector class.

#ifndef BITSET_H_
#define BITSET_H_

#include "mtl/Vec.h"
#include <limits>
namespace Monosat {

//=================================================================================================
#define BITSET_ELEMENT_SIZE (sizeof(uint64_t))
class Bitset {
    vec<uint64_t>  buf;
    int sz;

public:
    Bitset():sz(0)  {}
    Bitset(int size):sz(size)  {buf.growTo(size/BITSET_ELEMENT_SIZE + 1);}
    Bitset(int size, bool default_value):sz(size)  {
    	buf.growTo(size/BITSET_ELEMENT_SIZE + 1);
    	if(default_value){
    		for(int i = 0;i<buf.size();i++){
    			buf[i]= std::numeric_limits<uint64_t>::max();
    		}
    	}
    }

    void memset(bool to){
    	if(to){
			for(int i = 0;i<buf.size();i++){
				buf[i]= std::numeric_limits<uint64_t>::max();
			}
		}else{
			for(int i = 0;i<buf.size();i++){
				buf[i]= 0;
			}
		}
    }

    void growTo(int size){
    	buf.growTo(size/BITSET_ELEMENT_SIZE + 1);
    	sz=size;
    }
    void clear (bool dealloc = false) { buf.clear(dealloc);sz=0; }
    void zero(){
    	for(int i = 0;i<buf.size();i++)
    		buf[i]=0;
    }
    int  size  () const { return sz; }
    void copyFrom(const Bitset & from){
    	from.buf.copyTo(buf);
    	sz=from.size();
    }
    void copyTo(Bitset & to) const{
       	buf.copyTo(to.buf);
       	to.sz=sz;

       }
    inline const bool operator [] (int index) const  {
    	int i = index/BITSET_ELEMENT_SIZE ;
    	int rem = index % BITSET_ELEMENT_SIZE;
    	assert(i<buf.size());
    	return buf[i] & (1<<rem);
    }

    inline void set(int index){
    	assert(index<size());
    	int i = index/BITSET_ELEMENT_SIZE;
    	int r = index %BITSET_ELEMENT_SIZE;

    	buf[i]|=(1<<r);

    }
    inline void clear(int index){
    	assert(index<size());
    	int i = index/BITSET_ELEMENT_SIZE;
    	int r = index %BITSET_ELEMENT_SIZE;

    	buf[i]&= ~(1<<r);

    }
    inline void toggle(int index){
    	assert(index<size());
    	int i = index/BITSET_ELEMENT_SIZE;
    	int r = index %BITSET_ELEMENT_SIZE;

    	buf[i]^= (1<<r);

    }
    void Not(Bitset & out){
        	out.clear();
        	out.growTo(size());

        	for(int i = 0;i<buf.size();i++){
        		uint64_t a = buf[i];
        		out.buf[i]=~a;
        	}
        	out.sz=sz;
      }
    void And(const Bitset & with){

        	int max = size();
        	if(max>with.size()){
        		max= with.size();
        	}
        	int max_i=max/BITSET_ELEMENT_SIZE+1;
        	for(int i = 0;i<max_i;i++){

        		uint64_t b = with.buf[i];
        		buf[i]&=b;
        	}

        }

    void Or(const Bitset & with){

          	int max = size();
          	if(max>with.size()){
          		max= with.size();
          	}
          	int max_i=max/BITSET_ELEMENT_SIZE+1;
          	for(int i = 0;i<max_i;i++){

          		uint64_t b = with.buf[i];
          		buf[i]|=b;
          	}

          }
    void And(const Bitset & with, Bitset & out){
    	out.clear();
    	out.growTo(size());
    	int max = size();
    	if(max>with.size()){
    		max= with.size();
    	}
    	int max_i=max/BITSET_ELEMENT_SIZE+1;
    	for(int i = 0;i<max_i;i++){
    		uint64_t a = buf[i];
    		uint64_t b = with.buf[i];
    		out.buf[i]=a&b;
    	}
    	for(int i = max_i;i<buf.size();i++){
    		uint64_t a = buf[i];
    		out.buf[i]=a;
    	}
    	out.sz=sz;
    }

    void Or(const Bitset & with, Bitset & out){
        	out.clear();
        	out.growTo(size());
        	int max = size();
        	if(max>with.size()){
        		max= with.size();
        	}
        	int max_i=max/BITSET_ELEMENT_SIZE+1;
        	for(int i = 0;i<max_i;i++){
        		uint64_t a = buf[i];
        		uint64_t b = with.buf[i];
        		out.buf[i]=a|b;
        	}
        	for(int i = max_i;i<buf.size();i++){
        		uint64_t a = buf[i];
        		out.buf[i]=a;
        	}
        	out.sz=sz;
        }
};
//=================================================================================================
}

#endif /* BITSET_H_ */
