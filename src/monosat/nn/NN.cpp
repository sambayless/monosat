/*
 * NN.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: sam
 */

#include "monosat/nn/NN.h"
#include "caffe/util/math_functions.hpp"
#include <cstdio>
#include <fstream>
#include "monosat/mtl/Rnd.h"
#include <limits>

using namespace caffe;  // NOLINT(build/namespaces)
using std::string;

 NeuralNetwork::NeuralNetwork(const std::string & model_file,const std::string & trained_file,double rnd_seed):rnd_seed(rnd_seed){

		#ifdef CPU_ONLY
		Caffe::set_mode(Caffe::CPU);
		#else
		Caffe::set_mode(Caffe::GPU);
		#endif
		//following google deepdream, it looks like caffe doesn't currently have a way to programatically access
		//the force backward parameter
		string patched_model_file = std::tmpnam(nullptr);//fix this later...
		std::ifstream  src(model_file, std::ios::binary);
		std::ofstream  dst(patched_model_file,   std::ios::binary);

		dst << src.rdbuf();
		dst<< "\nforce_backward: true\n";
		dst.close();
		src.close();

		net =new Net<float>(patched_model_file, TEST);//not sure if this should be TEST or TRAIN
		std::remove(patched_model_file.c_str());

		net->CopyTrainedLayersFrom(trained_file);


		int n_outs=  net->num_outputs();
		int n_ins = net->num_inputs();
		assert(n_ins==1);
		assert(n_outs==1);

		 input_layer = net->input_blobs()[0];
		input_shape = input_layer->shape();
		num_channels = input_layer->channels();
		input_geometry = cv::Size(input_layer->width(), input_layer->height());

		 output_layer = net->output_blobs()[0];
		 output_shape = output_layer->shape();

		input_lower.resize(nInputs(),0.0f);
		input_upper.resize(nInputs(),1.0f);

		output_lower.resize(nOutputs(),0.0f);
		output_upper.resize(nOutputs(),1.0f);

		input_hit_lower.resize(nInputs(),false);
		input_hit_upper.resize(nInputs(),false);

		output_hit_lower.resize(nOutputs(),false);
		output_hit_upper.resize(nOutputs(),false);
	}



 void  NeuralNetwork::forward(){

 		  net->ForwardPrefilled();

 		  /* Copy the output layer to a std::vector */
/*
 		  const float* begin = output_layer->cpu_data();
 		  const float* end = begin + output_layer->channels();
 		  output =  std::vector<float>(begin, end);*/
 		 // return std::vector<float>(begin, end);

 	}
 void NeuralNetwork::backward(){
	 net->ClearParamDiffs();
	 net->Backward();

	/* src = net.blobs['data'] # input image is stored in Net's 'data' blob
	    dst = net.blobs[end]

	    ox, oy = np.random.randint(-jitter, jitter+1, 2)
	    src.data[0] = np.roll(np.roll(src.data[0], ox, -1), oy, -2) # apply jitter shift

	    net.forward(end=end)
	    objective(dst)  # specify the optimization objective
	    net.backward(start=end)
	    g = src.diff[0]
	    # apply normalized ascent step to the input image
	    src.data[:] += step_size/np.abs(g).mean() * g*/



 }
bool NeuralNetwork::solve(long limit){

	for (int i =0;i<nInputs();i++){
		input_hit_lower[i]=false;
		input_hit_upper[i]=false;
	}
	for (int i =0;i<nOutputs();i++){
		output_hit_lower[i]=false;
		output_hit_upper[i]=false;
	}

	long iterations=0;
	//first: put all inputs into the legal range (where this method will keep them)
	for(int i = 0;i<nOutputs();i++){
		if(output_lower[i]>output_upper[i]){
			return false;
		}
	}
	if(!initializeValues(true))
		return false;

#ifndef NDEBUG
		for(int i = 0;i<nInputs();i++){
			float cur_val = input_layer->cpu_data()[i] ;
			if(cur_val<input_lower[i]){
				assert(false);
			}else if(cur_val>input_upper[i]){
				assert(false);
			}
		}
#endif

	bool computed_once=false;

	bool all_inputs_fixed=true;
	for(int i = 0;i<nInputs();i++){

		if(fabs(input_upper[i]-input_lower[i])>0 ){
			all_inputs_fixed=false;
			break;
		}
	}
	bool all_out_of_range=false;
	while(true){
		forward();
		//Need to change the error signal here, to deal with ranges on the outputs.
		//Set the diff to the smallest value needed to get the output back within range.
		bool all_in_range=true;
		for(int i = 0;i<nOutputs();i++){
			float cur_val = output_layer->cpu_data()[i] ;
			float cur_diff = output_layer->cpu_diff()[i] ;
			if(cur_val<output_lower[i]){
				output_hit_lower[i]=true;
				float mid_point =output_lower[i]+( output_upper[i]-output_lower[i])/2.0;
				output_layer->mutable_cpu_diff()[i] =-(mid_point-cur_val); //or should the difference try to place the value in the mid range, rather than the end point, for faster convergence?
				assert(	output_layer->cpu_diff()[i]<0 );
				all_in_range=false;
			}else if(cur_val>output_upper[i]){
				output_hit_upper[i]=true;
				float mid_point =output_lower[i]+( output_upper[i]-output_lower[i])/2.0;
				output_layer->mutable_cpu_diff()[i] =-(mid_point- cur_val);
				assert(	output_layer->cpu_diff()[i]>0 );
				all_in_range=false;
			}else{
				output_layer->mutable_cpu_diff()[i]=0;//this output value is in the acceptable range.
			}
		}

		if(all_in_range){
			//all outputs are in range
#ifndef NDEBUG
			for(int i = 0;i<nOutputs();i++){
				float cur_val = output_layer->cpu_data()[i] ;
				if(cur_val<output_lower[i]- epsilon()){
					assert(false);
				}else if(cur_val>output_upper[i]+ epsilon()){
					assert(false);
				}
			}
			for(int i = 0;i<nInputs();i++){
				float cur_val = input_layer->cpu_data()[i] ;
				if(cur_val<input_lower[i]- epsilon()){
					assert(false);
				}else if(cur_val>input_upper[i]+ epsilon()){
					assert(false);
				}
			}
#endif
	/*		for(int i = 0;i<nOutputs();i++){
				float cur_val = output_layer->cpu_data()[i] ;
				printf("%d: %f, ",i,cur_val);
			}
			printf("\n");*/
			for(int i = 0;i<nOutputs();i++){
					float cur_val = output_layer->cpu_data()[i] ;
					printf(" %d: %f, ",i,cur_val);
				}
				printf("\n");
			int x =0;
			int y = 0;
			for(int i = 0;i<nInputs();i++){
				//printf("%d ",input_lower[i]>0.5 ? 1:0);
				printf("%.2f ",input_lower[i]);
				x++;
				if(x>=28){
					x=0;
					y++;
					printf("\n");
				}

			}
			printf("\n");
			return true;
		}else if (all_out_of_range){
			all_out_of_range=false;
			if(!initializeValues(false,InitPolicy::RANDOM))
				return false;
			continue;
		}

		backward();
#ifndef NDEBUG
		for(int i = 0;i<nInputs();i++){
			float cur_val = input_layer->cpu_data()[i] ;
			if(cur_val<input_lower[i]-epsilon()){
				assert(false);
			}else if(cur_val>input_upper[i]+epsilon()){
				assert(false);
			}
		}
#endif
		if(iterations++>limit){
			if(all_out_of_range){
				printf("All input errors out of range, giving up...");
			}
			for(int i = 0;i<nOutputs();i++){
				float cur_val = output_layer->cpu_data()[i] ;
				printf("%d: %f, ",i,cur_val);
			}
			printf("\n");
			break;//could also check for a fixed point and re-initialize values randomly in that case
		}
		if(all_inputs_fixed && computed_once)
			return false;
		all_out_of_range=true;
		//src.data[:] += step_size/np.abs(g).mean() * g*/

		//note:force diffs to stay within the allowed upper/lower bounds.

		/*const vector<float>& net_params_lr = this->net->params_lr();
		  float momentum = this->param_.momentum();
		  float local_rate = learning_rate * net_params_lr[param_id];*/
		  // Compute the update to history, then copy it to the parameter diff.

		  switch (Caffe::mode()) {
		  case Caffe::CPU: {
			  caffe_cpu_axpby(input_layer->count(), learning_rate,
		    		input_layer->cpu_diff(),1.0f, input_layer->mutable_cpu_diff());
			    //clip all the input values so they stay in range.
				for(int i = 0;i<nInputs();i++){
					float old_val = input_layer->cpu_data()[i] ;
					float diff = -input_layer->cpu_diff()[i];
					//float new_val = old_val+ diff;
					assert(input_upper[i]>=input_lower[i]);
					assert(old_val>=input_lower[i]-epsilon());
					assert(old_val<=input_upper[i]+epsilon());
					if(old_val+diff>input_upper[i]-epsilon()){
						input_hit_upper[i]=true;
						diff = input_upper[i]-old_val;
						assert(old_val+diff==input_upper[i]);
					}else if(old_val+diff<input_lower[i]+epsilon()){
						input_hit_lower[i]=true;
						diff = input_lower[i]-old_val;
						if(!floatEq(old_val+diff,input_lower[i])){
							bool a = floatEq(old_val+diff,input_lower[i]);
						}
						assert(floatEq(old_val+diff,input_lower[i]));
					}
					if(fabs(diff)>epsilon()){//or some small value..
						all_out_of_range=false;
					}

					input_layer->mutable_cpu_diff()[i]=-diff;
				}

		    break;
		  }
		  case Caffe::GPU: {
		#ifndef CPU_ONLY
			  caffe_gpu_axpby(input_layer->count(), learning_rate,
			  		    		input_layer->gpu_diff(),1.0f, input_layer->mutable_gpu_diff());
		#else
		    NO_GPU;
		#endif
		    break;
		  }
		  }
		  computed_once=true;
		input_layer->Update();//Note: The diff's are _negated_ before being added to data

		//important: No other layers or weights are updated in this process!
#ifndef NDEBUG
		for(int i = 0;i<nInputs();i++){
			float cur_val = input_layer->cpu_data()[i] ;
			if(cur_val<input_lower[i] - epsilon()){
				assert(false);
			}else if(cur_val>input_upper[i]+ epsilon()){
				assert(false);
			}
		}
#endif

	}
	return false;
}

float rnd_float(float lower, float upper,double & seed) {
	assert(lower<=upper);
	float rnd = Monosat::drand(seed);
    float diff = upper - lower;
    float r = rnd * diff;
    float v = lower + r;
    assert(v>=lower);
    assert(v<=upper);
    return v;
}

bool NeuralNetwork::initializeValues(bool only_out_of_bounds, InitPolicy init_policy){
	if (init_policy==InitPolicy::DEFAULT)
		init_policy = default_init_policy;
	for(int i = 0;i<nInputs();i++){
		float cur_val = input_layer->cpu_data()[i] ;
		float lower = input_lower[i];
		float upper = input_upper[i];
		if(lower>upper){
			return false;//this is unsat
		}else if(lower==upper){
			setModel(upper,i,true);
			continue;
		}
		if(only_out_of_bounds){
			switch(init_policy){
				case InitPolicy::NEAREST_BOUND:
					if(cur_val<lower){
						setModel(lower,i,true);
					}else if(cur_val>upper){
						setModel(upper,i,true);
					}
					break;
				case InitPolicy::FARTHEST_BOUND:
					if(cur_val<lower){
						setModel(upper,i,true);
					}else if(cur_val>upper){
						setModel(lower,i,true);
					}
					break;
				case InitPolicy::RANDOM:
					if(cur_val<lower || cur_val>upper){
						float n_val = rnd_float(lower,upper,rnd_seed);
						setModel(n_val,i,true);
					}
					break;
				case InitPolicy::MIDDLE:
					if(cur_val<lower || cur_val>upper){
						float v = lower + (upper - lower) / 2.0;
						setModel(v,i,true);
					}
					break;
			}
		}else{
			switch(init_policy){
				case InitPolicy::NEAREST_BOUND:
					setModel(lower,i,true);
					break;
				case InitPolicy::FARTHEST_BOUND:
					setModel(upper,i,true);
					break;
				case InitPolicy::RANDOM:

					setModel(rnd_float(lower,upper,rnd_seed),i,true);
					break;
				case InitPolicy::MIDDLE:

					setModel(lower + (upper - lower) / 2.0,i,true);
					break;
			}
		}
	}
	//correspondingly init the output values
	forward();
	return true;
}
