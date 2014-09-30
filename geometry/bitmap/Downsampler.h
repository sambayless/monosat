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
using namespace Monosat;

//Takes a high-dimensional generator and creates a version with a down-sampled and normalized parameter space
class Downsampler:public BitmapGenerator{
public:
	BitmapGenerator & g;
	int n_samples;
	int n_rnd_samples;
	vec<vec<int> > sampling_points;//the values at each parameter that are sampled
	void generate(Bitmap & bitmap){
		for(int p = 0;p<nParams();p++){
			int v = getParameter(p);
			assert(v<sampling_points[p].size());
			assert(v>=0);

			int sampling_point = sampling_points[p][v];
			g.setParameter(p,sampling_point);
		}
		g.generate(bitmap);
	}

	Downsampler(BitmapGenerator & _g, double & seed, bool only_downsample_seed=false):g(_g),n_samples(8), n_rnd_samples(8){
		if(g.hasRandomSeed()){
			int max = g.getMaximum(0);
			int samples = std::min(max,n_rnd_samples);
			addParameter(samples-1);
			sampling_points.push();
			if(samples==max+1){
				for(int i = 0;i<samples;i++){
					sampling_points[0].push(i);
				}
			}else{
				for(int i = 0;i<samples;i++){
					int sample = irand(seed,max+1);
					if(!sampling_points[0].contains(sample)){
						sampling_points[0].push(sample);
					}
				}
			}
		}

		for(int p = g.hasRandomSeed()? 1:0;p<g.nParams();p++){
			int max = g.getMaximum(p);
			int samples = std::min(max,n_samples);
			if(only_downsample_seed){
				samples= max+1;
			}
			sampling_points.push();
			if(g.isFixed(p)){
				addParameter(0);
				samples = 1;
				sampling_points[p].push(g.getParameter(p));
			}else{
				addParameter(samples-1);
				if(samples==max+1){
					for(int i = 0;i<samples;i++){
						sampling_points[p].push(i);
					}
				}else{
					for(int i = 0;i<samples;i++){
						int sample = irand(seed,max+1);
						if(!sampling_points[p].contains(sample)){
							sampling_points[p].push(sample);
						}
					}
				}
			}
		}
	}

};



#endif /* DOWNSAMPLER_H_ */
