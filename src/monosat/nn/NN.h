/*
 * NN.h
 *
 *  Created on: Nov 7, 2015
 *      Author: sam
 */

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <limits>
#define CPU_ONLY
#include <caffe/caffe.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>


class NeuralNetwork{
public:
	enum class InitPolicy{
		DEFAULT,NEAREST_BOUND,FARTHEST_BOUND,RANDOM, MIDDLE
	};
	InitPolicy default_init_policy = InitPolicy::NEAREST_BOUND;
	private:
	double rnd_seed;
	caffe::Net<float> * net=nullptr;
	cv::Size input_geometry;
	caffe::Blob<float>* input_layer;
	caffe::vector<int> input_shape;
	int num_channels=-1;
	caffe::Blob<float>* output_layer;
	caffe::vector<int> output_shape;

	float learning_rate=0.01;
	cv::Mat input;//Current input layer values

	std::vector<float> output;//Current output values
	std::vector<float> output_lower;
	std::vector<float> output_upper;
	std::vector<float> input_lower;
	std::vector<float> input_upper;

	//cv::Mat & output;//Current output layer values


	void WrapInputLayer(std::vector<cv::Mat>* input_channels);
	void Preprocess(const cv::Mat& img, std::vector<cv::Mat>* input_channels);
	//Feed forward the network from the inputs, updating the outputs
	void forward();

	void backward();//Must be called after 'forward'. Compute error between actual outputs and expected outputs, and backpropagate to get deltas on the inputs of the network
public:
	std::vector<bool> input_hit_upper;
	std::vector<bool> input_hit_lower;
	std::vector<bool> output_hit_upper;
	std::vector<bool> output_hit_lower;
	NeuralNetwork(const std::string &model_file,const std::string& trained_file,double rnd_seed = 1);
	//Iteratively modify the input data to maximimze
	bool solve(long limit=std::numeric_limits<long>::max());

	int nInputs(){
		return input_layer->count();
	}
	int nOutputs(){
		return output_layer->count();
	}
	int getIndex(std::vector<int> & indices,bool input){
		if(input){
			return input_layer->offset(indices);
		}else{
			return output_layer->offset(indices);
		}
	}
/*	float getInput(int inputID){
		int row = inputID/input.cols;
		int col = inputID % input.cols;
		return input.at<float>(row,col);
	}*/
	float getOutput(int outputID){
		return output_layer->cpu_data()[outputID];
		//return output[outputID];
	}
	float getOutputDelta(int outputID){
		return output_layer->cpu_diff()[outputID];
		//return output[outputID];
	}
	float getInput(int inputID){
		return input_layer->cpu_data()[inputID];
		//return input_layer-> data_at(n, channel,row, col);
	}

	//Get the error delta at an input
	float getInputDelta(int inputID){

		return input_layer->cpu_diff()[inputID];
		//return input_layer-> diff_at(n, channel,row, col);
	}

	float getModel(int index, bool input){
		if(input){
			return getInput(index);
		}else{
			return getOutput(index);
		}
	}

	float getBound(int index, bool input, bool lower){
		if(input){
			if(lower){
				return input_lower[index];
			}else{
				return input_upper[index];
			}
		}else{
			if(lower){
				return output_lower[index];
			}else{
				return output_upper[index];
			}
		}
	}


	void setBound(float value,int index, bool input, bool lower){
		if(input){
			if(lower){
				input_lower[index]=value;
			}else{
				input_upper[index]=value;
			}
		}else{
			if(lower){
				output_lower[index]=value;
			}else{
				output_upper[index]=value;
			}
		}
	}
	void setModel(float value,int index, bool input){
		if(input){
			input_layer->mutable_cpu_data()[index]=value;
		}else{
			output_layer->mutable_cpu_data()[index]=value;
		}
	}

	bool modelInBounds(int index,bool input){

		float model = getModel(index,input);
		if(model< getBound(index,input,true)-epsilon()){
			return false;
		}
		if(model> getBound(index,input,false)+epsilon()){
			return false;
		}
		return true;
	}

	bool initializeValues(bool only_out_of_bounds, InitPolicy init_policy=InitPolicy::DEFAULT);
	float epsilon()const{
		return 0.001;
	}
	bool floatEq(float a, float b){
		float diff = a-b;
		if(diff<0){
			return diff>-epsilon();
		}else{
			return diff<epsilon();
		}
	}
};

