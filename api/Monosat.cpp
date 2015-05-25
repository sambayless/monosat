
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Solver.h"
#include "simp/SimpSolver.h"
#include "graph/GraphTheory.h"
#include "geometry/GeometryTheory.h"
#include "pb/PbTheory.h"
#include "bv/BVTheorySolver.h"
#include "core/SolverTypes.h"
#include "Monosat.h"
#include "mtl/Vec.h"
using namespace Monosat;
using namespace std;

//Select which algorithms to apply for graph and geometric theory solvers, by parsing command line arguments and defaults.
void selectAlgorithms(){
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
		fprintf(stderr, "Error: unknown max-flow/min-cut algorithm %s, aborting\n",
				((string) opt_maxflow_alg).c_str());
		exit(1);
	}

	componentsalg = ComponentsAlg::ALG_DISJOINT_SETS;

	if (!strcasecmp(opt_components_alg, "disjoint-sets")) {
		componentsalg = ComponentsAlg::ALG_DISJOINT_SETS;
	} else {
		fprintf(stderr, "Error: unknown connectivity algorithm %s, aborting\n",
				((string) opt_components_alg).c_str());
		exit(1);
	}

	cyclealg = CycleAlg::ALG_DFS_CYCLE;

	if (!strcasecmp(opt_cycle_alg, "dfs")) {
		cyclealg = CycleAlg::ALG_DFS_CYCLE;
	}else if (!strcasecmp(opt_cycle_alg, "pk")) {
		cyclealg = CycleAlg::ALG_PK_CYCLE;
	} else {
		fprintf(stderr, "Error: unknown cycle detection algorithm %s, aborting\n",
				((string) opt_cycle_alg).c_str());
		exit(1);
	}


	mstalg = MinSpanAlg::ALG_KRUSKAL;

	if (!strcasecmp(opt_mst_alg, "kruskal")) {
		mstalg = MinSpanAlg::ALG_KRUSKAL;
	} else if (!strcasecmp(opt_mst_alg, "prim")) {
		mstalg = MinSpanAlg::ALG_PRIM;
	} else if (!strcasecmp(opt_mst_alg, "spira-pan")) {
		mstalg = MinSpanAlg::ALG_SPIRA_PAN;
	} else {
		fprintf(stderr, "Error: unknown minimum spanning tree algorithm %s, aborting\n",
				((string) opt_mst_alg).c_str());
		exit(1);
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
		fprintf(stderr, "Error: unknown reachability algorithm %s, aborting\n", ((string) opt_reach_alg).c_str());
		exit(1);
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
		fprintf(stderr, "Error: unknown distance algorithm %s, aborting\n", ((string) opt_dist_alg).c_str());
		exit(1);
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
		fprintf(stderr, "Error: unknown undirected reachability algorithm %s, aborting\n",
				((string) opt_reach_alg).c_str());
		exit(1);
	}

	allpairsalg = AllPairsAlg::ALG_DIJKSTRA_ALLPAIRS;

	if (!strcasecmp(opt_allpairs_alg, "floyd-warshall")) {
		allpairsalg = AllPairsAlg::ALG_FLOYDWARSHALL;
	} else if (!strcasecmp(opt_allpairs_alg, "dijkstra")) {
		allpairsalg = AllPairsAlg::ALG_DIJKSTRA_ALLPAIRS;

	} else {
		fprintf(stderr, "Error: unknown allpairs reachability algorithm %s, aborting\n",
				((string) opt_allpairs_alg).c_str());
		exit(1);
	}

	undirected_allpairsalg = AllPairsConnectivityAlg::ALG_DIJKSTRA_ALLPAIRS;

	if (!strcasecmp(opt_undir_allpairs_alg, "floyd-warshall")) {
		undirected_allpairsalg = AllPairsConnectivityAlg::ALG_FLOYDWARSHALL;
	} else if (!strcasecmp(opt_undir_allpairs_alg, "dijkstra")) {
		undirected_allpairsalg = AllPairsConnectivityAlg::ALG_DIJKSTRA_ALLPAIRS;
	} else if (!strcasecmp(opt_undir_allpairs_alg, "thorup")) {
		undirected_allpairsalg = AllPairsConnectivityAlg::ALG_THORUP;
	} else {
		fprintf(stderr, "Error: unknown undirected allpairs reachability algorithm %s, aborting\n",
				((string) opt_allpairs_alg).c_str());
		exit(1);
	}
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

#if defined(__linux__)
	fpu_control_t oldcw, newcw;
	_FPU_GETCW(oldcw);
	newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
	_FPU_SETCW(newcw);
	if (opt_verb > 0)
		fprintf(stderr, "WARNING: for repeatability, setting FPU to use double precision\n");
#endif
	selectAlgorithms();
	  Monosat::SimpSolver * S = new Monosat::SimpSolver();

	  S->_external_data =(void*)new MonosatData();

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
bool solve(Monosat::SimpSolver * S){
		S->preprocess();//do this _even_ if sat based preprocessing is disabled! Some of the theory solvers depend on a preprocessing call being made!
		S->eliminate(true);
	  static Monosat::vec<Monosat::Lit> assume;
	  assert(assume.size()==0);
	  return S->solve(assume);
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
bool solveAssumption(Monosat::SimpSolver * S,int * assumptions, int n_assumptions){
	  static Monosat::vec<Monosat::Lit> assume;
		S->preprocess();//do this _even_ if sat based preprocessing is disabled! Some of the theory solvers depend on a preprocessing call being made!
		S->eliminate(true);
	  assume.clear();
	  for (int i = 0;i<n_assumptions;i++){
		  assume.push(toLit(assumptions[i]));
	  }
	  return S->solve(assume);
 }

 int newVar(Monosat::SimpSolver * S){
	  return S->newVar();
 }

 bool addClause(Monosat::SimpSolver * S,int * assumptions, int n_assumptions){
	  static vec<Lit> clause;
	  clause.clear();
	  for (int i = 0;i<n_assumptions;i++){
		  clause.push(toLit(assumptions[i]));
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

 void bv_addition( Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID1, int bvID2, int resultID){
	  bv->newAdditionBV(resultID,bvID1,bvID2);
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


 //model query
 //Returns 0 for true, 1 for false, 2 for unassigned.
 int getModel_Literal(Monosat::SimpSolver * S,int lit){
	 return toInt(S->value(toLit(lit)));
 }
 long getModel_BV(Monosat::SimpSolver * S, Monosat::BVTheorySolver<long> * bv, int bvID){
	 return bv->getUnderApprox(bvID);
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

 long getModel_MinimumSpanningTreeWeight(Monosat::SimpSolver * S,Monosat::GraphTheorySolver<long> *G,int spanning_tree_literal){
	 Lit l = toLit(spanning_tree_literal);
	 return G->getModel_MaximumFlowModel(S->getTheoryLit(l));
 }
