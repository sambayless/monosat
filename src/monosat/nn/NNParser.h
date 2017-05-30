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

#ifndef NN_PARSER_H_
#define NN_PARSER_H_

#include <stdio.h>

#include "monosat/utils/ParseUtils.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/nn/NNTheory.h"
#include "monosat/core/Config.h"
#include "monosat/core/Dimacs.h"
#include <set>
#include <string>
#include <sstream>
namespace Monosat {


//=================================================================================================
// GRAPH Parser:
template<class B, class Solver>
class NNParser: public Parser<B, Solver> {
	using Parser<B, Solver>::mapVar;
	using Parser<B, Solver>::mapBV;
	vec<NNTheory<long>*> nns;
	vec<Var> vars;
	vec<char>  tmp;

	struct Comparison{
		int nnID;
		Lit l;
		bool input;
		int index;
		float geq;
	};
	vec<Comparison> comps;
	struct BVNN{
		int nnID;
		bool input;
		std::vector<int> indices;
		int bvID;

		float min;
		float max;
	};
	vec<BVNN> bvnns;
public:
	NNParser():Parser<B, Solver> ("Neural Network") {
		
	}

	void readFilename(B & in, string & s){
		skipWhitespace(in);
		while(!isEof(in) && *in != '\n' && *in != ';'){
			s+= *in;
			++in;
		}
		size_t endpos = s.find_last_not_of(" \t");
		if( string::npos != endpos )
		{
		    s = s.substr( 0, endpos+1 );
		}

	}

	bool parseLine(B& in, Solver& S) {
		
		skipWhitespace(in);
		if (*in == EOF)
			return false;
		if (match(in,"neural_network")){
			//neural_network <id> prototype_filename ; trained_weights_filename
			int id = parseInt(in);
			int semi_pos=-1;
			string prototype_name;
			readFilename(in,prototype_name);
			++in;
			string weights_name;
			readFilename(in,weights_name);
			nns.growTo(id+1,nullptr);
			if(nns[id]){
				parse_errorf("Re-defined neural network %d", id);
			}
			if(prototype_name.size()==0 || weights_name.size()==0){
				parse_errorf("Missing network definition");
			}
			//Make this portable!
			if(prototype_name[0] != '/' && prototype_name[0]!='~'){
				prototype_name = string(opt_nn_path)+ "/" + prototype_name;
			}
			if(weights_name[0] != '/' && weights_name[0]!='~'){
				weights_name = string(opt_nn_path)+ "/" + weights_name;
			}
			nns[id]=new NNTheory<long>(&S);
			nns[id]->loadNeuralNetwork(prototype_name,weights_name);
			return true;
		}else if (match(in,"nn")){
			//nn <nnID> input <Lit>  %d >= %f
			int id = parseInt(in);
			if(id>=nns.size() || !nns[id]){
				parse_errorf("Undefined neural network %d", id);
			}
			skipWhitespace(in);
			if (match(in,"bv")){
				skipWhitespace(in);
				int bvID = parseInt(in);
				skipWhitespace(in);
				bool input = false;
				if(match(in,"input")){
					input=true;
				}else if(!match(in,"output")){
					parse_errorf("bad nn bv\n");
				}
				int dimensions = parseInt(in);
				std::vector<int> indices;
				for(int i =0;i<dimensions;i++){
					int index = parseInt(in);
					indices.push_back(index);
				}


				float min = parseFloat(in,tmp);
				float max = parseFloat(in,tmp);
				bvnns.push({id,input,indices,bvID,min,max});
				return true;
			}else{
				if(match(in,"input")){
					int parsed_lit = parseInt(in);
					if (parsed_lit == 0)
						parse_errorf("invalid dimacs literal (0)");
					Var var = abs(parsed_lit) - 1;
					var= mapVar(S,var);
					Lit l = mkLit(var,parsed_lit<0);
					if(sign(l)){
						parse_errorf("neural network literals must be positive");
					}
					int index = parseInt(in);
					skipWhitespace(in);
					if(match(in,">=")){
						float geq = parseFloat(in,tmp);
						comps.push();
						comps.last().nnID = id;
						comps.last().l = l;
						comps.last().input=true;
						comps.last().index=index;
						comps.last().geq=geq;
					}else{
						parse_errorf("Unsupported neural network operation %s", in);
					}
					return true;
				}else if(match(in,"output")){
					int parsed_lit = parseInt(in);
					if (parsed_lit == 0)
						parse_errorf("invalid dimacs literal (0)");
					Var var = abs(parsed_lit) - 1;
					var= mapVar(S,var);
					Lit l = mkLit(var,parsed_lit<0);
					if(sign(l)){
						parse_errorf("neural network literals must be positive");
					}
					int index = parseInt(in);
					skipWhitespace(in);
					if(match(in,">=")){
						float geq = parseFloat(in,tmp);
						comps.push();
						comps.last().nnID = id;
						comps.last().l = l;
						comps.last().input=false;
						comps.last().index=index;
						comps.last().geq=geq;
					}else{
						parse_errorf("Unsupported neural network operation %s", in);
					}
					return true;
				}
			}
		}

		return false;
	}
	
	void implementConstraints(Solver & S) {
/*		for(auto & c: comps){
			int nnID =c.nnID;
			nns[nnID]->addGEQ(c.index,c.input,c.geq,var(c.l));
		}*/
		for(auto & b:bvnns){
			int nnID =b.nnID;
			b.bvID= mapBV(S,b.bvID);
			int bvID = b.bvID;

			BVTheorySolver<long>* bvTheory =(BVTheorySolver<long>*) S.getBVTheory();

			if(!bvTheory->hasBV(bvID)){
				parse_errorf("Undefined bvID");
			}
			nns[nnID]->newNNBv(b.input,b.indices,bvID,b.min,b.max);


		}
	}

	
};

//=================================================================================================
}
;

#endif /* GRAPH_PARSER_H_ */
