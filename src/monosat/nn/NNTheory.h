/****************************************************************************************[Solver.h]
 The MIT License (MIT)

 Copyright (c) 2015, Sam Bayless

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

#ifndef NNTHEORY_H_
#define NNTHEORY_H_

#include "monosat/mtl/Vec.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/core/Theory.h"
#include "monosat/core/Solver.h"
#include "monosat/bv/BVTheorySolver.h"
#include "monosat/bv/BVTheory.h"
#include "monosat/nn/NN.h"
#include "monosat/mtl/Rnd.h"
#include <limits>
#include <cmath>
namespace Monosat {

template<typename Weight>
class NNTheory: public Theory,public BVTheory {
	Solver * S;
	BVTheorySolver<Weight>* bvTheory;
	double rnd_seed;
	int theory_index=-1;
	bool connectedToBV=false;
	int nnid=-1;
public:
	int getNNID()const{
		return nnid;
	}
	void setNNID(int nnid){
		this->nnid = nnid;
	}
	CRef assign_reason;


	vec<Lit> tmp_clause;

	double propagationtime=0;
	long stats_propagations=0;
	long stats_propagations_skipped=0;
	long stats_shrink_removed = 0;
	long stats_reasons = 0;
	long stats_conflicts = 0;
	bool needs_propagation=false;

	NeuralNetwork * nn=nullptr;

	struct NNVar{
		bool input;
		int index;
	};


	struct BVData{
		bool input;
		int index=-1;
		float min;
		float max;
	};
	vec<BVData> bvs;

	vec<int> input_bvs;
	vec<int> output_bvs;

	vec<NNVar> to_propagate;
	vec<bool> inputs_to_propagate;
	vec<bool> outputs_to_propagate;
	int inconsistent_index=-1;
	int inconsistent_level=-1;
	bool inconsistent_is_input=false;
	bool inconsistent_model=true;

	CRef output_geq;
	CRef output_leq;

	int n_exact_input_assignments=0;

	struct Assignment{
		bool refined_under;

		int bvID;
		//float diff; //in principle, could just store diffs here, but I'm worried about error building up...
		//Actually, if a _binary_ diff representation is used here, rather than a floating point diff, then this might work.
		//bool changed;
		float prev;
		float next;
		float prev_model;
	};
	vec<Assignment> trail;
	vec<int> trail_lim;

	vec<bool> fixed_input;
public:
	void setBVTheory(BVTheorySolver<Weight> * bv){
		if(bvTheory!=bv){
			bvTheory = bv;

			if(bv){

				bv->addTheory(this);
			}
		}
	}
	NNTheory(Solver * S) :
			S(S){


		rnd_seed=drand(S->getRandomSeed());

		S->addTheory(this);
		output_geq = S->newReasonMarker(this);
		output_leq = S->newReasonMarker(this);
		//assign_false_reason=S->newReasonMarker(this);
		//assign_true_reason=S->newReasonMarker(this);

		//Atoms are of the form (output or input) >= constant
	}
	~NNTheory() {
	}
	;

	void set_needs_prop(bool input, int index){
		if(input){
			if(!inputs_to_propagate[index]){
				inputs_to_propagate[index]=true;
				to_propagate.push({input,index});
			}
		}else{
			if(!outputs_to_propagate[index]){
				outputs_to_propagate[index]=true;
				to_propagate.push({input,index});
			}
		}
	}

	void newNNBv(bool input,std::vector<int> & indices,int bvID,float min,float max){
		assert(bvTheory);
		if(!connectedToBV){
			connectedToBV=true;
			bvTheory->addTheory(this);

		}
		bvTheory->setBitvectorTheory(bvID,getTheoryIndex());
		bvs.growTo(bvID+1);
		assert(bvs[bvID].index==-1);
		int index = nn->getIndex(indices,input);assert(index>=0);
		bvs[bvID].index=index;
		bvs[bvID].input=input;
		bvs[bvID].min=min;
		bvs[bvID].max=max;
		if(input){
			input_bvs[index]=bvID;
		}else{
			output_bvs[index]=bvID;
		}
		nn->setBound(max,index,input,false);
		nn->setBound(min,index,input,true);

	}
	inline int getTheoryIndexBV(){
		return theory_index;
	}
/*	inline double getBucketSize(int bvID){
		BVData & data = bvs[bvID];
		BitVector<Weight> bv = bvTheory->getBV(bvID);
		double min = (double) data.min;
		double max = (double) data.max;
		int width = bv.width();
		long maxval = (1L<<width);
		double bucket_size = (max-min)/((double)maxval);
		return bucket_size;
	}*/
	inline float getFloat(int bvID, Weight val){
		BVData & data = bvs[bvID];
		BitVector<Weight>  bv = bvTheory->getBV(bvID);
		double min = (double) data.min;
		double max = (double) data.max;
		int width = bv.width();
		if(val==0){
			return min;
		}
		Weight maxval = (1L<<width)-1;
		if (val==maxval){
			return max;
		}
		double d = (max-min)/((double)maxval);
		return (float)((val*d) + min);
	}

	inline Weight discretize(float value,int width, float min,float max, bool under){
		if(value==min){
			return 0;
		}
		Weight maxval = (1L<<width)-1;
		if (value==max){
			return maxval;
		}
		double d = (max-min)/((double)maxval);
		if(under){
			return (Weight) floor((value-min)/d);
		}else{
			return (Weight) ceil((value-min)/d);
		}
	}

/*	inline float getUnder(int bvID, Weight val){
		BVData & data = bvs[bvID];
		BitVector<Weight>  bv = bvTheory->getBV(bvID);
		double min = (double) data.min;
		double max = (double) data.max;
		int width = bv.width();


		double bucket_size = getBucketSize(bvID);
		Weight bucket = bv.getUnder();
		return (float)( bucket_size*bucket) + min;
	}
	inline float getOver(int bvID, Weight val){
		BVData & data = bvs[bvID];
		BitVector<Weight> bv = bvTheory->getBV(bvID);
		double min = (double) data.min;
		double max = (double) data.max;
		int width = bv.width();
		long maxval = (1L<<width)-1;
		double bucket_size = getBucketSize(bvID);
		Weight bucket = bv.getOver();
		if(bucket<maxval)
			return (float)( bucket_size*(bucket+1)+ min - std::numeric_limits<float>::epsilon());
		else{
			return max;
		}
	}*/
	Weight epsilon()const{
		return 0.001;
	}
	void enqueueBV(int bvID){
		while (S->decisionLevel() > decisionLevel()) {
				newDecisionLevel();
			}
		if(bvs[bvID].index>-1){
			if(inconsistent_index<0){


				float prev_bound=-1;

				BVData & data = bvs[bvID];
				BitVector<Weight>  bv = bvTheory->getBV(bvID);
				Weight u = bv.getUnder();
				Weight o = bv.getOver();
				float under = getFloat(bvID, bv.getUnder());
				float over = getFloat(bvID, bv.getOver());

				if(data.input && data.index==22){
					int a=1;
				}

				if(under>  nn->getBound(data.index,data.input,true)){
					if(!data.input && data.index!=6){
						int a=1;
					}
					bool caused_refinement=false;
					prev_bound = nn->getBound(data.index,data.input,true);
					assert(under>prev_bound);
					caused_refinement=true;

					nn->setBound(under,data.index,data.input,true);

					if(under > nn->getBound(data.index,data.input,false)+epsilon()){
						//this would make the upper and lower bounds inconsistent
						inconsistent_index=data.index;
						inconsistent_is_input =data.input;
						inconsistent_level = decisionLevel();
						S->needsPropagation(getTheoryIndex());
					}else{
						//diff = nn->getBound(data.index,data.input,sign(l)) - data.geq;
						//assert(diff>0);
						if(!nn->modelInBounds(data.index,data.input)){
							inconsistent_model=true;
							S->needsPropagation(getTheoryIndex());
						}
						if(data.input && ! fixed_input[data.index]){
						if(fabs(nn->getBound(data.index,data.input,false)- nn->getBound(data.index,data.input,true)) ==0 ){
							n_exact_input_assignments++;
							fixed_input[data.index]=true;
							assert(n_exact_input_assignments<=nn->nInputs());
						}
						}
						//set_needs_prop(data.input,data.index);
					}

					if(caused_refinement){
						trail.push({true,bvID, prev_bound, under,nn->getModel(data.index,data.input)});

					}
				}

				if (over<  nn->getBound(data.index,data.input,false)){
					if(!data.input)
						return;//refuse to set upper bound on the output of the nn, for now
					bool caused_refinement=false;
					prev_bound = nn->getBound(data.index,data.input,false);
					assert(over<prev_bound);
					caused_refinement=true;

					nn->setBound(over,data.index,data.input,false);

					if(over < nn->getBound(data.index,data.input,true)-epsilon()){
						//this would make the upper and lower bounds inconsistent
						inconsistent_index=data.index;
						inconsistent_is_input =data.input;
						inconsistent_level = decisionLevel();
						S->needsPropagation(getTheoryIndex());
					}else{
						//diff = nn->getBound(data.index,data.input,sign(l)) - data.geq;
						//assert(diff>0);
						if(!nn->modelInBounds(data.index,data.input)){
							inconsistent_model=true;
							S->needsPropagation(getTheoryIndex());
						}
						if(data.input && ! fixed_input[data.index]){
						if(fabs(nn->getBound(data.index,data.input,false)- nn->getBound(data.index,data.input,true)) ==0 ){
							n_exact_input_assignments++;
							fixed_input[data.index]=true;
							assert(n_exact_input_assignments<=nn->nInputs());
						}
						}
						//set_needs_prop(data.input,data.index);
					}

					if(caused_refinement){
						trail.push({false,bvID, prev_bound, over,nn->getModel(data.index,data.input)});

					}
				}




			}


		}
	};
	void backtrackBV(int bvID){



	}

	void loadNeuralNetwork(const std::string & model_file,const std::string & trained_file){
		assert(!nn);
		nn = new NeuralNetwork(model_file,trained_file,drand(S->random_seed));
		if(opt_nn_init==0){
			nn->default_init_policy = NeuralNetwork::InitPolicy::NEAREST_BOUND;
		}else if(opt_nn_init==1){
			nn->default_init_policy = NeuralNetwork::InitPolicy::FARTHEST_BOUND;
		}else if(opt_nn_init==2){
			nn->default_init_policy = NeuralNetwork::InitPolicy::RANDOM;
		}else if(opt_nn_init==3){
			nn->default_init_policy = NeuralNetwork::InitPolicy::MIDDLE;
		}
		input_bvs.growTo(nn->nInputs()+1,-1);
		output_bvs.growTo(nn->nOutputs()+1,-1);
		fixed_input.growTo(nn->nInputs(),false);
	}



	inline int getTheoryIndex() const {
		return theory_index;
	}
	inline void setTheoryIndex(int id) {
		theory_index = id;
	}
	inline void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	inline int decisionLevel() {
		return trail_lim.size();
	}

	inline void backtrackUntil(int untilLevel){

		if(untilLevel<decisionLevel()){
			for (int c = trail.size() - 1; c >= trail_lim[untilLevel]; c--) {
				Assignment & e = trail[c];
				int bvID = e.bvID;
				BVData & data = bvs[bvID];
				if(data.input  && fixed_input[data.index]){
					if(fabs(nn->getBound(data.index,data.input,false)- nn->getBound(data.index,data.input,true)) ==0 ){
						if(e.prev != nn->getBound(data.index,data.input,e.refined_under)){
							n_exact_input_assignments--;
							fixed_input[data.index]=false;
							assert(n_exact_input_assignments>=0);
						}
					}
				}
				nn->setBound(e.prev,data.index,data.input,e.refined_under);
				nn->setModel(e.prev_model,data.index,data.input);
			}
			//inconsistent_model=false;
			trail.shrink(trail.size() - trail_lim[untilLevel]);
			trail_lim.shrink(trail_lim.size() - untilLevel);
			if(decisionLevel()<inconsistent_level){
				assert(inconsistent_index>-1);
				inconsistent_level = -1;
				inconsistent_index=-1;
				inconsistent_is_input=false;
			}
		}
	}

	inline void undecideTheory(Lit l){

	}

	void enqueueTheory(Lit l) {


	}
	vec<Lit> tmp_conflict;
	vec<int> rnd_inputs;
	void buildNaiveClause(vec<Lit> & conflict,bool input){
		if(!input){
			for(int i =0;i<nn->nOutputs();i++){
				int bvID = output_bvs[i];
				bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);
				bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
			}
		}else{
			for(int i = 0;i<nn->nInputs();i++){
				int bvID = input_bvs[i];
				bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);
				bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
			}
		}
		//bvTheory->toSolver(conflict);
	}

	void buildBoundsClause(vec<Lit> & conflict, bool input){
		if(!input){
			for(int i =0;i<nn->nOutputs();i++){

				if(!nn->modelInBounds(i,false)){
					int bvID = output_bvs[i];
					float model = nn->getModel(i,false);
					if(model>nn->getBound(i,false,false)){
						//the upper bound has to be relaxed
						bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);
					}else{
						assert(model<nn->getBound(i,false,true));
						//the lower bound has to be relaxed
						bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
					}
				}
			}
		}else{
			for(int i = 0;i<nn->nInputs();i++){
				if(!nn->modelInBounds(i,true)){
					int bvID = input_bvs[i];
					float model = nn->getModel(i,true);
					if(model>nn->getBound(i,true,false)){
						//the upper bound has to be relaxed
						bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);
					}else{
						assert(model<nn->getBound(i,true,true));
						//the lower bound has to be relaxed
						bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
					}
				}
			}
		}
		//bvTheory->toSolver(conflict);
	}

	void buildBackPropClause(vec<Lit> & conflict,bool input){
		if(!input){
			//build an approximate learnt clause from the current linear approximation of the neural network
			//(aka, from the back propagation error signal)

			//for any variables where the error signal would put the variable out of its bounds, add that bound to the conflict.
			for(int i =0;i<nn->nOutputs();i++){
				if(!nn->modelInBounds(i,false)){
					int bvID = output_bvs[i];
					float model = nn->getModel(i,false);
					if(model>nn->getBound(i,false,false)){
						//the upper bound has to be relaxed
						//Lit r = output_upper_bound_reason[i];
						bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);

					}else{
						assert(model<nn->getBound(i,false,true));
						//the lower bound has to be relaxed
						bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
					}
				}
			}
		}else{

			for(int i = 0;i<nn->nInputs();i++){
				int bvID = input_bvs[i];
				float model = nn->getModel(i,true);
				float error = -nn->getInputDelta(i);
				if(model+error > nn->getBound(i,true,false)){
					bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);

				}else if(model+error < nn->getBound(i,true,true)){
					//the lower bound has to be relaxed
					bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
				}

			}

		}
		//bvTheory->toSolver(conflict);
	}
	std::vector<std::pair<float,int>> clause;

	void buildBackPropClause_max(vec<Lit> & conflict,int max_lits,bool input){
		assert(max_lits>0);

		   clause.clear();
			if(!input){
				//build an approximate learnt clause from the current linear approximation of the neural network
				//(aka, from the back propagation error signal)
				throw std::runtime_error("Unimplemented learning function");
				//for any variables where the error signal would put the variable out of its bounds, add that bound to the conflict.
				assert(false);
			}else{
				std::vector<char> changes;
				changes.clear();
				for(int i = 0;i<nn->nInputs();i++){
					int bvID = input_bvs[i];
					float model = nn->getModel(i,true);
					float error = -nn->getInputDelta(i);
					float diff = 0;
					changes.push_back(nn->getBound(i,true,true) >0.5 ? '1':'0');
					if(model+error > nn->getBound(i,true,false)){
						diff = model+error - nn->getBound(i,true,false);
						if(bvTheory->getOverApprox(bvID)== bvTheory->getOverApprox(bvID,true)){
							diff=0;
						}
						//bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);

					}else if(model+error < nn->getBound(i,true,true)){
						//the lower bound has to be relaxed
						diff = nn->getBound(i,true,true) - model+error ;
						if(bvTheory->getUnderApprox(bvID)== bvTheory->getUnderApprox(bvID,true)){
							diff=0;
						}
						//bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
					}
					assert(diff>=0);
					if(diff>0){
						clause.push_back({diff,i});
					}
				}
				std::sort(clause.begin(),clause.end(),[]( const pair<float,int> &a, const pair<float,int> &b ){ return a.first > b.first; } );
				for(int j = 0;j<std::min((int)clause.size(),max_lits);j++){
					int i = clause[j].second;
					int bvID = input_bvs[i];
					float model = nn->getModel(i,true);
					float error = -nn->getInputDelta(i);
					float diff = 0;
					if(model+error > nn->getBound(i,true,false)){
						changes[i]='+';
						bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);
					}else if(model+error < nn->getBound(i,true,true)){
						//the lower bound has to be relaxed
						changes[i]='-';
						bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
					}
				}
				int x =0;
				int y = 0;
				if(opt_verb>1){
				for(int i = 0;i<nn->nInputs();i++){
					printf("%c ",changes[i]);
					x++;
					if(x>=28){
						x=0;
						y++;
						printf("\n");
					}

				}
				}
			}
			//bvTheory->toSolver(conflict);
		}
	void buildBoundHitClause(vec<Lit> & conflict, bool input){
		if(!input){
			//build an approximate learnt clause from the current linear approximation of the neural network
			//(aka, from the back propagation error signal)

			//for any variables where the error signal would put the variable out of its bounds, add that bound to the conflict.
 			for(int i =0;i<nn->nOutputs();i++){
 				int bvID = output_bvs[i];
				bool hit_lower = nn->output_hit_lower[i];
				bool hit_upper = nn->output_hit_upper[i];
				if(hit_upper){
					bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);
				}
				if(hit_lower){
					//the lower bound has to be relaxed
					bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
				}

			}
		}else{
			int x =0;
			int y = 0;

			for(int i = 0;i<nn->nInputs();i++){
				int bvID = input_bvs[i];
				float model = nn->getModel(i,true);
				bool hit_lower = nn->input_hit_lower[i];
				bool hit_upper = nn->input_hit_upper[i];
				char change = nn->getBound(i,true,true) >0.5 ? '1':'0';
				if(hit_upper){

					bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),conflict);
					if(bvTheory->getOverApprox(bvID)!= bvTheory->getOverApprox(bvID,true)){
						change='+';
					}
				}
				if(hit_lower){
					//the lower bound has to be relaxed
					bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),conflict);
					if(bvTheory->getUnderApprox(bvID)!= bvTheory->getUnderApprox(bvID,true)){
						change='-';
					}
				}
				if(opt_verb>1){
				printf("%c ",change);

				x++;
				if(x>=28){
					x=0;
					y++;
					printf("\n");
				}
				}
			}
		}
		//bvTheory->toSolver(conflict);
	}

	void dbg_sync(){
//#ifndef NDEBUG
/*		for(int i = 0;i<nn->nInputs();i++){
			float lower = nn->getBound(i,true,true);
			float upper = nn->getBound(i,true,false);

		}*/

		int count_fixed = 0;
		for(int i = 0;i<nn->nInputs();i++){
			if(nn->getBound(i,true,true)== nn->getBound(i,true,false)){
				count_fixed++;
				if(!fixed_input[i]){
					throw std::runtime_error("Bad fixed input status!");
				}
			}else if(fixed_input[i]){
				throw std::runtime_error("Bad fixed input status (disabled)!");
			}
		}
		if(count_fixed!= n_exact_input_assignments){
			throw std::runtime_error("Bad fixed input count!");
		}
/*
		for(int bvID = 0;bvID<bvs.size();bvID++){
			if(bvs[bvID].index>-1){
					if(inconsistent_index<0){


						float prev_bound=-1;

						BVData & data = bvs[bvID];
						int i = data.index;
						BitVector<Weight>  bv = bvTheory->getBV(bvID);
						Weight u = bv.getUnder();
						Weight o = bv.getOver();
						float under = getFloat(bvID, bv.getUnder());
						float over = getFloat(bvID, bv.getOver());
						float lower = nn->getBound(i,data.input,true);
						float upper = nn->getBound(i,data.input,false);
						assert(under==lower);
						assert(over==upper);
						if(upper!=over){
							fprintf(stderr,"not in sync\n");exit(1);
						}
						if(under!=lower){
							fprintf(stderr,"not in sync\n");exit(1);
						}
					}
			}
		}*/
//#endif
	}

	//If nn only has linear and relu neurons, then it should be possible to do completely precise clause learning here
	//but there is no guarantee that the clauses will be small enough to be useful...
	bool propagateTheory(vec<Lit> & conflict)override {
		bool do_solve = false;
		bool do_enqueue=false;
		if(opt_nn_prop_during){
			do_solve |=(n_exact_input_assignments==nn->nInputs());
			do_enqueue |= n_exact_input_assignments==nn->nInputs();
		}
		return propagateTheory(conflict,do_solve,do_enqueue);
	}
	bool propagateTheory(vec<Lit> & conflict, bool do_solve, bool do_model_enqueue) {
		bool any_changed_model=false;
		int nin = nn->nInputs();
		static int iter=0;
		static int outer_it=0;
		dbg_sync();
		//technically, should also to basic propagation of bounds
		if(inconsistent_model && do_solve){
			++outer_it;
			if(outer_it==27){
				int a=1;
			}
			if(nn->getBound(4,false,false)<1.0f){
				int a=1;
			}
			if(opt_verb>1)
				printf("Fully solving nn %d...\n",outer_it);
			if(!nn->solve(opt_nn_backprop_limit)){
				++iter;
				if(opt_verb>1)
					printf("Gave up solving nn %d...\n",iter);
				if(iter==81){
					int a=1;
				}

				//if the network failed to converge to a satisfying model after a reasonable time, then treat this
				//as a conflict (even though a satisfying solution might in fact exist)

				//More clause learning options
				//1) Unconstrained deep dreaming
				//Deep dream without input constraints, producing a local optimum
				//Then either use the diff between that local optimum and the bounds to learn a clause,
				//Or the list of all bounds ever crossed while finding that optimum to learn the clause

				//2) local neighbourhood
				//Relax all bounds slightly, temporarily, during deep dreaming.
				//Then record which of those _relaxed_ bounds was hit during deep dreaming.
				//This gives the benefits of deep dreamining for clause learning, even if fully precise bounds are used.

				//also, option to union in all the decisions in the solver (or all the relevant decisions??)
				//In order to make the learnt clause apply only to the local part of the search space...

				if(opt_nn_learn_out==0){
					buildBoundHitClause(conflict,false);
				}else if (opt_nn_learn_out==1){
					if(opt_nn_learn_max<=0)
						buildBackPropClause(conflict,false);
					else
						buildBackPropClause_max(conflict,opt_nn_learn_max,false);
				}else if (opt_nn_learn_out==2){
					buildBoundsClause(conflict,false);
				}else if (opt_nn_learn_out==3){
					buildNaiveClause(conflict,false);
				}

				int output_lits =conflict.size();
				tmp_conflict.clear();
				if(opt_nn_learn_in==0){
					buildBoundHitClause(tmp_conflict,true);
				}else if (opt_nn_learn_in==1){
					if(opt_nn_learn_max<=0)
						buildBackPropClause(tmp_conflict,true);
					else
						buildBackPropClause_max(tmp_conflict,opt_nn_learn_max,true);
				}else if (opt_nn_learn_in==2){
					buildBoundsClause(tmp_conflict,true);
				}else if (opt_nn_learn_in==3){
					buildNaiveClause(tmp_conflict,true);
				}
				int input_lits = tmp_conflict.size();
				for (Lit l:tmp_conflict)
					conflict.push(l);
				tmp_conflict.clear();
				if(opt_nn_weaken_learnt_rnd>0){
					if(rnd_inputs.size()==0){
						for(int i =0;i<this->input_bvs.size();i++){
							if(input_bvs[i]>-1){
								rnd_inputs.push(i);
							}
						}
					}
					randomShuffle(rnd_seed, rnd_inputs);
					for(int i = 0;i<std::min((int)opt_nn_weaken_learnt_rnd,rnd_inputs.size());i++){
						int index = rnd_inputs[i];
						int bvID = input_bvs[index];
						if(bvID>=0){
							bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),tmp_conflict);
							bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),tmp_conflict);
						}
					}
					//bvTheory->toSolver(tmp_conflict);
					for (Lit l:tmp_conflict)
						conflict.push(l);
					tmp_conflict.clear();
				}
				if(opt_verb>1)
					printf("NN learnt clause size %d (%d input, %d output)\n",conflict.size(), input_lits,output_lits);
				return false;
			}else{
				printf("Solved nn...\n");
				inconsistent_model=false;
			}
		}else{
			//current model in the neural network still satisfies the constraints
		}
		if(do_model_enqueue && ! inconsistent_model){
			//can propagate the outputs of the neural network, since the inputs are fully specified.
			for(int i =0;i<output_bvs.size();i++){
				int bvID = output_bvs[i];
				if(bvID>=0){
					BVData & data = bvs[bvID];
					BitVector<Weight>  bv = bvTheory->getBV(bvID);
					assert(!data.input);
					float model = nn->getModel(data.index,data.input);
					assert(model>=nn->getBound(data.index,data.input,true));
					assert(model<=nn->getBound(data.index,data.input,false));
					Weight under = discretize(model,bv.width(), data.min,data.max,true);
					float under_m = getFloat(bvID,under);
					assert(under_m<= model);
					if(under_m>model)
						throw std::runtime_error("Inconsistent nn model propagation");

					if(under>bv.getUnder()){
						Lit l = bvTheory->toSolver(bvTheory->newComparison(Comparison::geq,bvID,under,var_Undef,false));
						if(S->value(l)==l_Undef){
							S->enqueue(l,output_geq);
						}else if (S->value(l)==l_False){
							buildReason(l,conflict,output_geq);
							return false;
						}
					}

					Weight over = discretize(model,bv.width(), data.min,data.max,false);
					float over_m = getFloat(bvID,over);
					assert(over_m>= model);
					if(over_m<model)
						throw std::runtime_error("Inconsistent nn model propagation");
					if(over<bv.getOver()){
						Lit l =  bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,over,var_Undef,false));
						if(S->value(l)==l_Undef){
							S->enqueue(l,output_leq);
						}else if (S->value(l)==l_False){
							buildReason(l,conflict,output_leq);
							return false;
						}
					}
				}
			}
		}
		return true;
	}
	void printStats(int detailLevel) {

	}
	void preprocess(){
		if(nn){
			nn->initializeValues(false);
			inconsistent_model=true;//after re-initializing the model values, mark the model as inconsistent.
		}
	}
	inline bool solveTheory(vec<Lit> & conflict){
		//inconsistent_model=true;
		//printf("Solving nn...\n");
		return propagateTheory(conflict,true,true);
	}
	inline void buildReason(Lit p, vec<Lit> & reason, CRef reason_marker){
		stats_reasons++;
		reason.push(bvTheory->toSolver(p));
		if(reason_marker==output_geq || reason_marker==output_leq){
			//just build naive clause for now... should really do proper clause learning here
			for(int i = 0;i<nn->nInputs();i++){
				int bvID = input_bvs[i];
				float model = nn->getModel(i,true);
				bool hit_lower = nn->input_hit_lower[i];
				bool hit_upper = nn->input_hit_upper[i];
				bvTheory->buildComparisonReason(Comparison::leq,bvID,bvTheory->getOverApprox(bvID),reason);
				//the lower bound has to be relaxed
				bvTheory->buildComparisonReason(Comparison::geq,bvID,bvTheory->getUnderApprox(bvID),reason);
			}
		}

	}
	bool check_solved() {
		for(int i =0;i<nn->nOutputs();i++){
			if(!nn->modelInBounds(i,false)){
				return false;
			}
		}

		for(int i = 0;i<nn->nInputs();i++){
			if(!nn->modelInBounds(i,true)){
				return false;
			}
		}
		return true;
	}
private:

	
};

}
;

#endif /* AMOTheory_H_ */
