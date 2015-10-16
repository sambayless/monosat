
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Solver.h"
#include "core/Config.h"
#include "simp/SimpSolver.h"
#include "graph/GraphTheory.h"
#include "geometry/GeometryTheory.h"
#include "pb/PbTheory.h"
#include "bv/BVTheorySolver.h"
#include "amo/AMOTheory.h"
#include "core/SolverTypes.h"
#include "Monosat.h"
#include "mtl/Vec.h"
#include "core/Dimacs.h"
#include "bv/BVParser.h"
#include "graph/GraphParser.h"
#include "utils/ParseUtils.h"
#include "amo/AMOParser.h"
#include <stdexcept>
#include <cstdarg>
using namespace Monosat;
using namespace std;



//Supporting function for throwing parse errors
inline void api_errorf(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char buf[1000];
    vsnprintf(buf, sizeof buf,fmt, args);
    va_end(args);
    fprintf(stderr, buf);
    fflush(stderr);
    throw std::runtime_error(buf);

}

//Select which algorithms to apply for graph and geometric theory solvers, by parsing command line arguments and defaults.
void _selectAlgorithms(){
	mincutalg = MinCutAlg::ALG_EDMONSKARP;

	if (!strcasecmp(opt_maxflow_alg, "edmondskarp-adj")) {
		mincutalg = MinCutAlg::ALG_EDKARP_ADJ;
	} else if (!strcasecmp(opt_maxflow_alg, "edmondskarp")) {
		mincutalg = MinCutAlg::ALG_EDMONSKARP;
	} else if (!strcasecmp(opt_maxflow_alg, "edmondskarp-dynamic")) {
		mincutalg = MinCutAlg::ALG_EDKARP_DYN;
	} else if (!strcasecmp(opt_maxflow_alg, "dinics")) {
		//Dinitz is also commonly spelled 'dinics' or 'Dinits', so accept those too...
		mincutalg = MinCutAlg::ALG_DINITZ;
	} else if (!strcasecmp(opt_maxflow_alg, "dinics-linkcut")) {
		mincutalg = MinCutAlg::ALG_DINITZ_LINKCUT;
	} else if (!strcasecmp(opt_maxflow_alg, "dinitz")) {
		mincutalg = MinCutAlg::ALG_DINITZ;
	} else if (!strcasecmp(opt_maxflow_alg, "dinitz-linkcut")) {
		mincutalg = MinCutAlg::ALG_DINITZ_LINKCUT;
	} else if (!strcasecmp(opt_maxflow_alg, "dinits")) {
		//Dinitz is also commonly spelled 'dinics' or 'Dinits', so accept those too...
		mincutalg = MinCutAlg::ALG_DINITZ;
	} else if (!strcasecmp(opt_maxflow_alg, "dinits-linkcut")) {
		mincutalg = MinCutAlg::ALG_DINITZ_LINKCUT;
	} else if (!strcasecmp(opt_maxflow_alg, "kohli-torr")) {
		mincutalg = MinCutAlg::ALG_KOHLI_TORR;
	} else {
		api_errorf( "Error: unknown max-flow/min-cut algorithm %s, aborting\n",
				((string) opt_maxflow_alg).c_str());

	}

	componentsalg = ComponentsAlg::ALG_DISJOINT_SETS;

	if (!strcasecmp(opt_components_alg, "disjoint-sets")) {
		componentsalg = ComponentsAlg::ALG_DISJOINT_SETS;
	} else {
		api_errorf(  "Error: unknown connectivity algorithm %s, aborting\n",
				((string) opt_components_alg).c_str());

	}

	cyclealg = CycleAlg::ALG_DFS_CYCLE;

	if (!strcasecmp(opt_cycle_alg, "dfs")) {
		cyclealg = CycleAlg::ALG_DFS_CYCLE;
	}else if (!strcasecmp(opt_cycle_alg, "pk")) {
		cyclealg = CycleAlg::ALG_PK_CYCLE;
	} else {
		api_errorf( "Error: unknown cycle detection algorithm %s, aborting\n",
				((string) opt_cycle_alg).c_str());

	}


	mstalg = MinSpanAlg::ALG_KRUSKAL;

	if (!strcasecmp(opt_mst_alg, "kruskal")) {
		mstalg = MinSpanAlg::ALG_KRUSKAL;
	} else if (!strcasecmp(opt_mst_alg, "prim")) {
		mstalg = MinSpanAlg::ALG_PRIM;
	} else if (!strcasecmp(opt_mst_alg, "spira-pan")) {
		mstalg = MinSpanAlg::ALG_SPIRA_PAN;
	} else {
		api_errorf( "Error: unknown minimum spanning tree algorithm %s, aborting\n",
				((string) opt_mst_alg).c_str());

	}

	reachalg = ReachAlg::ALG_BFS;

	if (!strcasecmp(opt_reach_alg, "dijkstra")) {
		reachalg = ReachAlg::ALG_DIJKSTRA;

	} else if (!strcasecmp(opt_reach_alg, "bfs")) {
		reachalg = ReachAlg::ALG_BFS;

	} else if (!strcasecmp(opt_reach_alg, "dfs")) {
		reachalg = ReachAlg::ALG_DFS;
	} else if (!strcasecmp(opt_reach_alg, "cnf")) {
		reachalg = ReachAlg::ALG_SAT;
	} else if (!strcasecmp(opt_reach_alg, "ramal-reps")) {
		reachalg = ReachAlg::ALG_RAMAL_REPS;
	} else {
		api_errorf( "Error: unknown reachability algorithm %s, aborting\n", ((string) opt_reach_alg).c_str());

	}

	distalg = DistAlg::ALG_DISTANCE;

	if (!strcasecmp(opt_dist_alg, "dijkstra")) {
		distalg = DistAlg::ALG_DIJKSTRA;

	} else if (!strcasecmp(opt_dist_alg, "bfs")) {
		distalg = DistAlg::ALG_DISTANCE;

	} else if (!strcasecmp(opt_dist_alg, "cnf")) {
		distalg = DistAlg::ALG_SAT;
	} else if (!strcasecmp(opt_dist_alg, "ramal-reps")) {
		distalg = DistAlg::ALG_RAMAL_REPS;
	} else {
		api_errorf(  "Error: unknown distance algorithm %s, aborting\n", ((string) opt_dist_alg).c_str());

	}

	undirectedalg = ConnectivityAlg::ALG_BFS;

	if (!strcasecmp(opt_con_alg, "dijkstra")) {
		undirectedalg = ConnectivityAlg::ALG_DIJKSTRA;

	} else if (!strcasecmp(opt_con_alg, "bfs")) {
		undirectedalg = ConnectivityAlg::ALG_BFS;

	} else if (!strcasecmp(opt_con_alg, "dfs")) {
		undirectedalg = ConnectivityAlg::ALG_DFS;
	} else if (!strcasecmp(opt_con_alg, "cnf")) {
		undirectedalg = ConnectivityAlg::ALG_SAT;
	} else if (!strcasecmp(opt_con_alg, "thorup")) {
		undirectedalg = ConnectivityAlg::ALG_THORUP;
	}/*else if (!strcasecmp(opt_con_alg,"ramal-reps")){
	 undirectedalg = ConnectivityAlg::ALG_RAMAL_REPS;
	 } */else {
		 api_errorf(  "Error: unknown undirected reachability algorithm %s, aborting\n",
				((string) opt_reach_alg).c_str());

	}

	allpairsalg = AllPairsAlg::ALG_DIJKSTRA_ALLPAIRS;

	if (!strcasecmp(opt_allpairs_alg, "floyd-warshall")) {
		allpairsalg = AllPairsAlg::ALG_FLOYDWARSHALL;
	} else if (!strcasecmp(opt_allpairs_alg, "dijkstra")) {
		allpairsalg = AllPairsAlg::ALG_DIJKSTRA_ALLPAIRS;

	} else {
		api_errorf(  "Error: unknown allpairs reachability algorithm %s, aborting\n",
				((string) opt_allpairs_alg).c_str());

	}

	undirected_allpairsalg = AllPairsConnectivityAlg::ALG_DIJKSTRA_ALLPAIRS;

	if (!strcasecmp(opt_undir_allpairs_alg, "floyd-warshall")) {
		undirected_allpairsalg = AllPairsConnectivityAlg::ALG_FLOYDWARSHALL;
	} else if (!strcasecmp(opt_undir_allpairs_alg, "dijkstra")) {
		undirected_allpairsalg = AllPairsConnectivityAlg::ALG_DIJKSTRA_ALLPAIRS;
	} else if (!strcasecmp(opt_undir_allpairs_alg, "thorup")) {
		undirected_allpairsalg = AllPairsConnectivityAlg::ALG_THORUP;
	} else {
		api_errorf(  "Error: unknown undirected allpairs reachability algorithm %s, aborting\n",
				((string) opt_allpairs_alg).c_str());

	}
}
void printStats(SimpSolver* solver) {
	double cpu_time = cpuTime();
    double mem_used = memUsedPeak(); // not available in osx

	solver->printStats(3);
	if (mem_used != 0)
		printf("Memory used           : %.2f MB\n", mem_used);
	printf("CPU time              : %g s\n", cpu_time);
}

struct MonosatData{
	Monosat::BVTheorySolver<long> * bv_theory=nullptr;
	vec< Monosat::GraphTheorySolver<long> *> graphs;
};
void * newSolver(int argc, char**argv){
	parseOptions(argc, argv, true);
	if (opt_adaptive_conflict_mincut == 1) {
		opt_conflict_min_cut = true;
		opt_conflict_min_cut_maxflow = true;
	}
	Monosat::opt_record=strlen(opt_record_file)>0;
	if(strlen(opt_debug_learnt_clauses)>0){
		opt_write_learnt_clauses=fopen(opt_debug_learnt_clauses,"w");
	}else{
		opt_write_learnt_clauses=nullptr;
	}
	_selectAlgorithms();
	  Monosat::SimpSolver * S = new Monosat::SimpSolver();

	  S->_external_data =(void*)new MonosatData();
	  if(!opt_pre){
		  S->eliminate(true);//disable preprocessing.
	  }
    return S ;
}
void deleteSolver (Monosat::SimpSolver * S)
{
	  if(S->_external_data){
		  delete((MonosatData*)S->_external_data);
		  S->_external_data=nullptr;
	  }
     delete (S);
}

void readGNF(Monosat::SimpSolver * S, const char  * filename){
	bool precise = true;

	gzFile in = gzopen(filename, "rb");
	if (in == nullptr)
		throw std::runtime_error("ERROR! Could not open file");


	Dimacs<StreamBuffer, SimpSolver> parser;
	BVParser<char *, SimpSolver> bvParser;
	parser.addParser(&bvParser);

	SymbolParser<char*,SimpSolver> symbolParser;
	parser.addParser(&symbolParser);

	GraphParser<char *, SimpSolver> graphParser(precise,bvParser.theory);
	parser.addParser(&graphParser);

	AMOParser<char *, SimpSolver> amo;
	parser.addParser(&amo);

	StreamBuffer strm(in);
	vec<int> assumps;
	bool ran_last_solve=false;
	while(parser.parse(strm, *S)){
		assumps.clear();
		for(Lit l:parser.assumptions){
			assumps.push(toInt(l));
		}
		solveAssumptions_MinBVs(S,&assumps[0],assumps.size(),&parser.bv_minimize[0],parser.bv_minimize.size());
		if(*strm==EOF){
			ran_last_solve=true;
		}
	}
	assert(*strm==EOF);
	if(!ran_last_solve){
		for(Lit l:parser.assumptions){
			assumps.push(toInt(l));
		}
		solveAssumptions_MinBVs(S,&assumps[0],assumps.size(),&parser.bv_minimize[0],parser.bv_minimize.size());
	}


	gzclose(in);

}

void * newGraph(Monosat::SimpSolver * S){
	  MonosatData * d = (MonosatData*) S->_external_data;
	  Monosat::GraphTheorySolver<long> *graph = new Monosat::GraphTheorySolver<long>(S);
	  S->addTheory(graph);
	  d->graphs.push(graph);
	  if( d->bv_theory){
		  graph->setBVTheory(d->bv_theory);
	  }

	  return graph;
}
void backtrack(Monosat::SimpSolver * S){
	S->cancelUntil(0);
}
void * initBVTheory(Monosat::SimpSolver * S){
	MonosatData * d = (MonosatData*) S->_external_data;
	if(d->bv_theory)
		return d->bv_theory;

	  Monosat::BVTheorySolver<long> * bv = new Monosat::BVTheorySolver<long>(S);

	  d->bv_theory=bv;
	  for (auto graph:d->graphs)
		  graph->setBVTheory(bv);

	  return bv;
}
bool solve(Monosat::SimpSolver * S){
	return solveAssumptions(S,nullptr,0);
  }

bool solveAssumptions(Monosat::SimpSolver * S,int * assumptions, int n_assumptions){
	return solveAssumptions_MinBVs(S,assumptions,n_assumptions,nullptr,0);
 }

bool solveAssumptions_MinBVs(Monosat::SimpSolver * S,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs){
	lbool ret = toLbool(solveAssumptionsLimited_MinBVs(S,-1,assumptions, n_assumptions, minimize_bvs, n_minimize_bvs));
	if(ret==l_True){
		return true;
	}else if (ret==l_False){
		return false;
	}else{
		throw std::runtime_error("Failed to solve!");
	}
}

int solveLimited(Monosat::SimpSolver * S,int time_cutoff){
	return solveAssumptionsLimited(S,time_cutoff,nullptr,0);
  }

int solveAssumptionsLimited(Monosat::SimpSolver * S,int time_cutoff,int * assumptions, int n_assumptions){
	return solveAssumptionsLimited_MinBVs(S,time_cutoff,assumptions,n_assumptions,nullptr,0);
 }

int solveAssumptionsLimited_MinBVs(Monosat::SimpSolver * S,int time_cutoff,int * assumptions, int n_assumptions, int * minimize_bvs, int n_minimize_bvs){
	using namespace Monosat;
	bool hit_cutoff=false;
	S->cancelUntil(0);

	  S->preprocess();//do this _even_ if sat based preprocessing is disabled! Some of the theory solvers depend on a preprocessing call being made!

	  vec<Monosat::Lit> assume;
	  for (int i = 0;i<n_assumptions;i++){
		   Lit l =toLit( assumptions[i]);
		   if (var(l)>=S->nVars()){
			   api_errorf("Assumption ltieral %d is not allocated",dimacs(l));
		   }
		  assume.push(l);
		  //S->setFrozen(v,true); //this is done in the solve() call
	  }

/*	  if (opt_pre){
		S->eliminate(false);//should this really be set to disable future preprocessing here?
	 }*/

	  vec<int> bvs;//bit vectors to minimize
	  for (int i = 0;i<n_minimize_bvs;i++){
		  int bvID = minimize_bvs[i];
		  if(!S->getBVTheory()){
			  api_errorf("No bitvector theory created (call initBVTheory())!");
		  }
		  if(! ((Monosat::BVTheorySolver<long> *) S->getBVTheory())->hasBV(bvID)){
			  api_errorf("Minimization bitvector %d is not allocated",bvID);
		  }
		  bvs.push(bvID);
	  }

	  bool r;
	  if (time_cutoff<0){
		  r= S->solve(assume,opt_pre,!opt_pre);
	  }else{
		  lbool res = S->solveLimited(assume,opt_pre,!opt_pre);
		  if (res==l_True){
			  r=true;
		  }else if (res==l_False){
			  r=false;
		  }else{
			  r=false;
			  hit_cutoff=true;

		  }
	  }
	  for(Lit l:assume){
	  		if(S->value(l)!=l_True){
	  			throw std::runtime_error("Model is inconsistent with assumptions!");
	  		}
	  	}
	  if(r && bvs.size()){
		  Monosat::BVTheorySolver<long> * bvTheory = (Monosat::BVTheorySolver<long> *) S->getBVTheory();
		  vec<long> min_values;
		  min_values.growTo(bvs.size());
		  long n_solves = 1;
		  for (int i = 0;i<bvs.size();i++){

			  int bvID = bvs[i];
			  long value = bvTheory->getUnderApprox(bvID);
			  if(opt_verb>=1){
				  printf("Minimizing bv%d (%d of %d)\n",bvID,i+1,bvs.size());
			  }

			  if(opt_verb>=1){
				  printf("Min bv%d = %ld",bvID,value);
			  }
			  min_values[i]=value;
			  long last_decision_value=value;
			  // int bvID,const Weight & to, Var outerVar = var_Undef, bool decidable=true
			  Lit last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,value,var_Undef,false));
			  while(value>0){
				  Lit decision_lit = bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,value-1,var_Undef,false));
				  assume.push(decision_lit);
				  n_solves++;

				  bool r;
				  if (time_cutoff<0){
					  r= S->solve(assume,false,false);
				  }else{
					  lbool res = S->solveLimited(assume,false,false);
					  if (res==l_True){
						  r=true;
					  }else if (res==l_False){
						  r=false;
					  }else{
						  r= false;
						  hit_cutoff=true;
					  }
				  }

				  if (r){
					  last_decision_lit=decision_lit;
					  last_decision_value=value-1;
					  if(S->value(decision_lit)!=l_True){
						  throw std::runtime_error("Error in optimization (comparison not enforced)");
					  }
					  for(Lit l:assume){
					  		if(S->value(l)!=l_True){
					  			throw std::runtime_error("Model is inconsistent with assumptions!");
					  		}
					  	}
					  long value2 = bvTheory->getUnderApprox(bvID);
					  if(value2>=value){
							throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");

					  }
					  value=value2;
					  min_values[i]=value;
					  assume.pop();
					  if(opt_verb>=1){
						  printf("\rMin bv%d = %ld",bvID,value);
					  }
				  }else{
					  assume.pop();
					  assert(min_values[i]==value);
					  assert(min_values[i]<=last_decision_value);
					  if(value<last_decision_value){
						  //if that last decrease in value was by more than 1
						  last_decision_lit =  bvTheory->toSolver(bvTheory->newComparison(Comparison::leq,bvID,value,var_Undef,false));
						  last_decision_value=value;
					  }
					  assume.push(last_decision_lit);
					  //this can be improved...
					  if (time_cutoff<0){
						  r= S->solve(assume,false,false);
					  }else{
						  lbool res = S->solveLimited(assume,false,false);
						  if (res==l_True){
							  r=true;
						  }else if (res==l_False){
							  r=false;
						  }else{
							  r= false;
							  hit_cutoff=true;
						  }
					  }
					  if(!r){
						  throw std::runtime_error("Error in optimization (instance has become unsat)");
					  }
					  for(Lit l:assume){
					  		if(S->value(l)!=l_True){
					  			throw std::runtime_error("Model is inconsistent with assumptions!");
					  		}
					  	}
					  if(min_values[i]< bvTheory->getUnderApprox(bvID)){
						  throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
					  }
					  break;
				  }
			  }
			  if(opt_verb>=1){
				  printf("\rMin bv%d = %ld\n",bvID, min_values[i]);
			  }

		  }
		  assert(r);

		  if(opt_verb>0){
			  printf("Minimum values found (after %ld calls) : ",n_solves);
			  for(int i = 0;i<min_values.size();i++){
				  int bvID = bvs[i];
				  printf("bv%d=%ld,",bvID,min_values[i]);
			  }
			  printf("\n");
		  }
		  if(opt_check_solution){
			  for(int i = 0;i<bvs.size();i++){
				  int bvID = bvs[i];
				  long min_value = min_values[i];
				  long model_val = bvTheory->getUnderApprox(bvID);
				  if(min_value!=model_val){
					  throw std::runtime_error("Error in optimization (minimum values are inconsistent with model)");
				  }
			  }
		  }
	  }
	if (opt_verb >= 1) {
		printStats(S);

	}
	if (hit_cutoff)
		 return toInt(l_Undef);
	return r?  toInt(l_True): toInt(l_False);
}

 int newVar(Monosat::SimpSolver * S){
	  return S->newVar();
 }

 bool disallowLiteralSimplification(Monosat::SimpSolver * S, int lit){
	 if(S->isEliminated(var(toLit(lit)))){
		 fprintf(stderr,"Warning: Literal %d has already been eliminated by the pre-processor\n", dimacs(toLit(lit)));
		 return false;
	 }else
		 S->setFrozen(var(toLit(lit)),true);
	 return true;
 }
 void setDecisionVar(Monosat::SimpSolver * S,int var,bool decidable){
	 S->setDecisionVar(var,decidable);
 }
 void setDecisionPriority(Monosat::SimpSolver * S,int var, int priority){
	 S->setDecisionPriority(var,priority);
 }
 bool isDecisionVar(Monosat::SimpSolver * S,int var){
	 return S->isDecisionVar(var);
 }
 int getDecisionPriority(Monosat::SimpSolver * S,int var){
	 return S->getDecisionPriority(var);
 }

 void setDecisionPolarity(Monosat::SimpSolver * S,Var v, bool b){
	  S->setPolarity(v,b);
 }

 bool getDecisionPolarity(Monosat::SimpSolver * S,Var v){
	return S->getPolarity(v);
 }


 int nVars(Monosat::SimpSolver * S){
	 return S->nVars();
 }

 int nClauses(Monosat::SimpSolver * S){
	 return S->nClauses();
 }
 int nBitvectors(Monosat::SimpSolver * S,Monosat::BVTheorySolver<long> * bv){
	 return bv->nBitvectors();
 }

 int true_lit(Monosat::SimpSolver * S){
	 return toInt(S->True());
 }
 bool addClause(Monosat::SimpSolver * S,int * lits, int n_lits){
	  static vec<Lit> clause;
	  clause.clear();
	  for (int i = 0;i<n_lits;i++){
		  clause.push(toLit(lits[i]));
	  }
	  return S->addClause(clause);
 }
 bool addUnitClause(Monosat::SimpSolver * S,int lit){
	  return S->addClause(toLit(lit));
 }
 bool addBinaryClause(Monosat::SimpSolver * S,int lit1, int lit2){
	  return S->addClause(toLit(lit1),toLit(lit2));
 }
 bool addTertiaryClause(Monosat::SimpSolver * S,int lit1, int lit2, int lit3){
	  return S->addClause(toLit(lit1),toLit(lit2),toLit(lit3));
 }

 //theory interface for bitvectors

 int newBitvector_const(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvWidth, long constval){
	 return bv->newBitvector(-1,bvWidth,constval).getID();
 }
 int newBitvector(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int * bits, int n_bits){
	  static vec<Var> lits;
	  lits.clear();
	  for (int i = 0;i<n_bits;i++){
		  lits.push(Var(bits[i]));
	  }
	  int bvID = bv->nBitvectors();
	  bv->newBitvector(bvID,lits);
	  return bvID;
 }

 int newBVComparison_const_lt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparison(Monosat::Comparison::lt,bvID,weight,v);
	  return toInt(l);
 }
 int newBVComparison_bv_lt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparisonBV(Monosat::Comparison::lt,bvID,compareID,v);
	  return toInt(l);
 }
 int newBVComparison_const_leq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparison(Monosat::Comparison::leq,bvID,weight,v);
	  return toInt(l);
 }
 int newBVComparison_bv_leq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparisonBV(Monosat::Comparison::leq,bvID,compareID,v);
	  return toInt(l);
 }

 int newBVComparison_const_gt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparison(Monosat::Comparison::gt,bvID,weight,v);
	  return toInt(l);
 }
 int newBVComparison_bv_gt(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparisonBV(Monosat::Comparison::gt,bvID,compareID,v);
	  return toInt(l);
 }
 int newBVComparison_const_geq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparison(Monosat::Comparison::geq,bvID,weight,v);
	  return toInt(l);
 }
 int newBVComparison_bv_geq(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, int compareID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  bv->newComparisonBV(Monosat::Comparison::geq,bvID,compareID,v);
	  return toInt(l);
 }
 void bv_min(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int n_args, int* args,int resultID){
	 vec<int> m_args;
	 for (int i = 0;i<n_args;i++)
		 m_args.push(args[i]);
	 bv->newMinBV(resultID, m_args);
 }
 void bv_max(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int n_args, int* args, int resultID){
	 vec<int> m_args;
	 for (int i = 0;i<n_args;i++)
		 m_args.push(args[i]);
	 bv->newMaxBV(resultID, m_args);
 }

 void bv_addition( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID1, int bvID2, int resultID){
	  bv->newAdditionBV(resultID,bvID1,bvID2);
 }
 void bv_subtraction( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID1, int bvID2, int resultID){
 	  bv->newSubtractionBV(resultID,bvID1,bvID2);
  }

void bv_ite( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int condition_lit,int bvThenID, int bvElseID, int bvResultID){
	Lit l = toLit(condition_lit);
	bv->newConditionalBV(l,bvThenID,bvElseID,bvResultID);
}

void bv_not(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int a,  int out){
	//return bv->bitwiseAnd(bv->getBV(a),bv->getBV(b)).getID();
	bv->bitwiseNot(bv->getBV(a),bv->getBV(out));
}

void bv_and(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int a, int b, int out){
	//return bv->bitwiseAnd(bv->getBV(a),bv->getBV(b)).getID();
	bv->bitwiseAnd(bv->getBV(a),bv->getBV(b),bv->getBV(out));
}
void bv_nand( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int a, int b, int out){
	//return bv->bitwiseNand(bv->getBV(a),bv->getBV(b)).getID();
	bv->bitwiseNand(bv->getBV(a),bv->getBV(b),bv->getBV(out));
}
void bv_or( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int a, int b, int out){
	//return bv->bitwiseOr(bv->getBV(a),bv->getBV(b)).getID();
	bv->bitwiseOr(bv->getBV(a),bv->getBV(b),bv->getBV(out));
}
void bv_nor( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int a, int b, int out){
	//return bv->bitwiseNor(bv->getBV(a),bv->getBV(b)).getID();
	bv->bitwiseNor(bv->getBV(a),bv->getBV(b),bv->getBV(out));
}
void bv_xor( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int a, int b, int out){
	//return bv->bitwiseXor(bv->getBV(a),bv->getBV(b)).getID();
	bv->bitwiseXor(bv->getBV(a),bv->getBV(b),bv->getBV(out));
}
void bv_xnor( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int a, int b, int out){
	//return bv->bitwiseXnor(bv->getBV(a),bv->getBV(b)).getID();
	bv->bitwiseXnor(bv->getBV(a),bv->getBV(b),bv->getBV(out));
}

void bv_concat( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int aID, int bID, int resultID){
	bv->concat(bv->getBV(aID), bv->getBV(bID),bv->getBV(resultID));
}

void bv_slice( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv,int aID, int lower, int upper, int resultID){
	bv->slice(bv->getBV(aID),lower,upper,bv->getBV(resultID));
}

 //simple at-most-one constraint: asserts that at most one of the set of variables (NOT LITERALS) may be true.
 //for small numbers of variables, consider using a direct CNF encoding instead
 void at_most_one(Monosat::SimpSolver * S, int * vars, int n_vars){
	 if(n_vars>1){
		 AMOTheory* amo = new  AMOTheory(S);
		 for(int i = 0;i<n_vars;i++){
			 Var v = vars[i];
			 amo->addVar(v);
		 }
	 }
 }

 //theory interface for graphs

 int newNode(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G){
	  return G->newNode();
 }
 int newEdge(Monosat::SimpSolver * S, Monosat::GraphTheorySolver<long> *G,int from,int  to,  long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->newEdge( from,  to, v,  weight );
	  return toInt(l);
 }
 int newEdge_bv(Monosat::SimpSolver * S, Monosat::GraphTheorySolver<long> *G,int from,int  to, int bvID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->newEdgeBV( from,  to, v,  bvID );
	  return toInt(l);
 }

 int reaches(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->reaches(from, to, v);
	  G->implementConstraints();
	  return toInt(l);
 }
 int shortestPathUnweighted_lt_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int steps){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->reaches(from, to, v,steps-1);
	  G->implementConstraints();
	  return toInt(l);
 }
 int shortestPathUnweighted_leq_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int steps){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->reaches(from, to, v,steps);
	  G->implementConstraints();
	  return toInt(l);
 }
 int shortestPath_lt_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, long dist){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->distance(from, to, v,dist, false);
	  G->implementConstraints();
	  return toInt(l);
 }
 int shortestPath_leq_const(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, long dist){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->distance(from, to, v,dist, true);
	  G->implementConstraints();
	  return toInt(l);
 }
 int shortestPath_lt_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int bvID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->distanceBV(from,to, v, bvID,false);
	  G->implementConstraints();
	  return toInt(l);
 }
 int shortestPath_leq_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int from, int to, int bvID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->distanceBV(from,to, v, bvID,true);
	  G->implementConstraints();
	  return toInt(l);
 }
 int maximumFlow_geq(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->maxflow(source, sink, v, weight,true);
	  G->implementConstraints();
	  return toInt(l);
 }
 int maximumFlow_gt(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->maxflow(source, sink, v, weight,false);
	  G->implementConstraints();
	  return toInt(l);
 }
 int maximumFlow_geq_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, int bvID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->maxflowBV(source, sink, v, bvID,true);
	  G->implementConstraints();
	  return toInt(l);
 }
 int maximumFlow_gt_bv(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, int bvID){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->maxflowBV(source, sink, v, bvID,false);
	  G->implementConstraints();
	  return toInt(l);
 }
 int minimumSpanningTree_leq(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->minimumSpanningTree(v, weight,true);
	  G->implementConstraints();
	  return toInt(l);
 }
 int minimumSpanningTree_lt(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int source, int sink, long weight){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->minimumSpanningTree(v, weight,false);
	  G->implementConstraints();
	  return toInt(l);
 }
 int acyclic_undirected(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->acyclic(v,false);
	  G->implementConstraints();
	  return toInt(l);
 }
 int acyclic_directed(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G){
	  Var v = newVar(S);
	  Lit l =mkLit(v);
	  G->acyclic(v,true);
	  G->implementConstraints();
	  return toInt(l);
 }

 void newEdgeSet(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int * edges, int n_edges){
	  static vec<int> edge_set;
	  edge_set.clear();
	  for (int i = 0;i<n_edges;i++){
		  Var outer=var(toLit(edges[i]));
		  if(outer>=S->nVars()){
			  api_errorf("Bad edge set variable %d",outer+1);
		  }
		  if(!S->hasTheory(outer)){
			  api_errorf("Bad edge set variable %d",outer+1);
		  }
		  if(S->getTheoryID(outer)!=G->getTheoryIndex()){
			  api_errorf("Wrong graph (%d) for variable %d",G->getTheoryIndex(),outer+1);
		  }
		  Var v = S->getTheoryVar(outer);
		  if(!G->isEdgeVar(v)){
			  api_errorf("Variable %d is not an edge variable",outer+1);
		  }
		  edge_set.push(G->getEdgeID(v));
	  }
	  G->newEdgeSet(edge_set);
 }

 //model query
 //Returns 0 for true, 1 for false, 2 for unassigned.
 int getModel_Literal(Monosat::SimpSolver * S,int lit){
	 Lit l = toLit(lit);
	 if (var(l)>=S->model.size())
		 return toInt(l_Undef);
	 lbool val = S->model[var(l)];
	 assert(val==l_True || val==l_False || val==l_Undef);
	 if (sign(l)){
		 if (val==l_True)
			 val=l_False;
		 else if (val==l_False){
			 val=l_True;
		 }
	 }
	 return toInt(val);//toInt(S->value(toLit(lit)));
 }
 long getModel_BV(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID, bool getMaximumValue){
	 if(getMaximumValue){
		 return bv->getOverApprox(bvID);
	 }else{
		 return bv->getUnderApprox(bvID);
	 }

 }
 //graph queries:
 long getModel_MaxFlow(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int maxflow_literal){
	 Lit l = toLit(maxflow_literal);
	 return G->getModel_MaximumFlow(S->getTheoryLit(l));
 }
 long getModel_EdgeFlow(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int maxflow_literal, int edge_assignment_literal){
	 Lit l = toLit(maxflow_literal);
	 Lit e = toLit(edge_assignment_literal);
	 return G->getModel_MaximumFlow_EdgeFlow(S->getTheoryLit(l),S->getTheoryLit(e));
 }
 long getModel_AcyclicEdgeFlow(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int maxflow_literal, int edge_assignment_literal){
	 Lit l = toLit(maxflow_literal);
	 Lit e = toLit(edge_assignment_literal);
	 return G->getModel_MaximumFlow_AcyclicEdgeFlow(S->getTheoryLit(l),S->getTheoryLit(e));
 }

 long getModel_MinimumSpanningTreeWeight(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int spanning_tree_literal){
	 Lit l = toLit(spanning_tree_literal);
	 return G->getModel_MinimumSpanningWeight(S->getTheoryLit(l));
 }
/* //Get the length of a valid path (from a reachability or shortest path constraint)
 long getModel_PathLength(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int reach_or_shortest_path_lit){
	 Lit l = toLit(reach_or_shortest_path_lit);
	 return G->getModel_PathLength(S->getTheoryLit(l));
 }
 //Get a valid path (from a reachability or shortest path constraint)
 //store_path must point to an array of ints of sufficient length to store the path (the path length can be optained by a call to getModel_PathLength)
 void getModel_Path(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int reach_or_shortest_path_lit, int * store_path){
	 Lit l = toLit(reach_or_shortest_path_lit);
	 std::vector<int> path;
	 G->getModel_Path(S->getTheoryLit(l),path);
	 for (int i = 0;i<path->size();i++){
		 store_path[i]=path[i];
	 }
	 //Returns the number of nodes in the path length for this reachability or shortest path literal (1+number of edges)
	   int getModel_PathLength(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int reach_or_shortest_path_lit){

	   }

		bool getModel_Path(Lit theoryLit,  int * store_path){

		 }
		 //Get a valid path, in terms of edges, (from a reachability or shortest path constraint)
		 //store_path must point to an array of ints of sufficient length to store the path (the path length can be optained by a call to getModel_PathLength)
		//Or, return false if there is no such path
		bool getModel_PathByEdgeLit(Lit theoryLit, int * store_path){

		}
 }

*/
