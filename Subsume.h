/*
 The MIT License (MIT)

 Copyright (c) 2013-2017, Sam Bayless

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
 */

#ifndef SUBSUME_H_
#define SUBSUME_H_

#include "core/SolverTypes.h"
#include "mtl/Vec.h"
#include "mtl/Sort.h"
namespace Minisat{
class Subsume{

	vec<int> seen;
	vec<uint32_t> abstractions;
	vec<int> clause_indices;
	int nvars = 0;


	struct LessThanSize{
		vec<vec<Lit> > & clauses;
		bool operator () (const int x, const int y){

			return clauses[x].size()<clauses[y].size();
		}

		LessThanSize(vec<vec<Lit> > & _clauses):clauses(_clauses){}
	};

	void sort_by_size(vec<int>& indices, vec<vec<Lit> > & clauses){
		sort(indices,LessThanSize(clauses));
	}

	public:

	inline bool subsumes(vec<Lit> & check, vec<Lit> & other,uint32_t abstract, uint32_t other_abstract) const
	{

	    if (other.size() < check.size() || (abstract & ~other_abstract) != 0)
	        return false;

	    bool subsumed=false;



	    for (unsigned i = 0; i < check.size(); i++) {
	        // search for c[i] or ~c[i]
	        for (unsigned j = 0; j < other.size(); j++)
	            if (check[i] == other[j])
	                goto ok;
	            else if (check[i] == ~(other[j])){
	            	return false;
	            }

	        // did not find it
	        return false;
	    ok:;
	    }

	    return true;
	}
	void setNumVars(int _nvars){
		nvars=_nvars;
	}

	void checkAll(vec<vec<Lit> > & attempt_to_subsume){
		abstractions.clear();
		clause_indices.clear();



		for(int i = 0;i<attempt_to_subsume.size();i++)
		{
			uint32_t a = calcAbstraction(attempt_to_subsume[i]);
			abstractions.push(a);
			clause_indices.push(i);
		}
		int num_subsumed=0;
		sort_by_size(clause_indices,attempt_to_subsume);
		for(int i = 0;i<clause_indices.size();i++){
			int clauseI = clause_indices[i];
			uint32_t a = abstractions[clauseI];
			int z = i+1;
			int j=i+1;
			for(;j<clause_indices.size();j++){
				int clauseJ = clause_indices[j];
				if(!subsumes(attempt_to_subsume[clauseI], attempt_to_subsume[clauseJ],a,abstractions[clauseJ])){
					clause_indices[z++]=clauseJ;
				}else{
					num_subsumed++;
				}
			}
			clause_indices.shrink(j-z);

		}
		if(num_subsumed>0){
			vec<vec<Lit> > tmp;

			for(int i = 0;i<clause_indices.size();i++){
				tmp.push();
				attempt_to_subsume[clause_indices[i]].copyTo(tmp.last());
			}

			for(int i = 0;i<tmp.size();i++){
				tmp[i].copyTo(attempt_to_subsume[i]);
			}
			attempt_to_subsume.shrink(attempt_to_subsume.size()-tmp.size());
		}


		printf("Subsumed %d clauses\n", num_subsumed);
	}

	void removeSubsumed(vec<Lit> & check, vec<vec<Lit> > & attempt_to_subsume, int start_index=0){


			seen.clear();
				seen.growTo(nvars);

					for(int j = 0;j<check.size();j++){
						Lit l = check[j];
						seen[toInt(l)]=true;
					}
					//uint32_t abstraction = calcAbstraction(check);
					//check if any interpolants in the super solver are subsumed:
					int j, h = start_index;
					int sz =attempt_to_subsume.size();
					for(j=start_index;j<attempt_to_subsume.size();j++){
						vec<Lit> & c =attempt_to_subsume[j];
						bool subsumed= false;



							if(check.size() && c.size()>check.size()){ //&&  ((abstraction & ~c.abstraction()) == 0)) {
								int count = 0;
								for(int k = 0;k<c.size();k++){// && (gen.size()- count > (c.size()-k) )
									if(seen[toInt(c[k])]){
										count++;
									}
								}

								assert(count<=check.size());
								subsumed= count==check.size();

							}



						if(subsumed){

						}else{
							//This is super inefficient!
							if(h!=j){
							   attempt_to_subsume[j].copyTo(attempt_to_subsume[h]);
							}
							h++;
						}
					}
					attempt_to_subsume.shrink(j-h);
	}
};
};

#endif /* SUBSUME_H_ */
