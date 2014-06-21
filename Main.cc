/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

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
#include <fstream>
#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include <signal.h>
#include <zlib.h>
#include <sstream>
#include <string>
#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "Aiger.h"
#include "graph/GraphParser.h"
#include "core/Dimacs.h"
#include "core/AssumptionParser.h"
#include "core/Solver.h"
#include "Aiger.h"
#include "core/Config.h"
#include <unistd.h>
#include <sys/time.h>
#include <algorithm>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include "simp/SimpSolver.h"
#include "pb/PbTheory.h"
#include "pb/PbParser.h"
#include "mtl/Map.h"
#include "geometry/GeometryTheory.h"
using namespace Minisat;
using namespace std;
//=================================================================================================


void printStats(Solver& solver)
{
		double cpu_time = cpuTime();
	    double mem_used = memUsedPeak();
    solver.printStats(3);
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
   printf("CPU time              : %g s\n", cpu_time);
}


static SimpSolver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt();
printf("\n"); printf("*** INTERRUPTED ***\n");
if (opt_verb > 0){
    printStats(*solver);
    printf("\n"); printf("*** INTERRUPTED ***\n"); }
fflush(stdout);
_exit(1);
}

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
/*    if (opt_verb > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }*/
    fflush(stdout);
    _exit(1); }


//=================================================================================================
// Main:
//from http://stackoverflow.com/a/994647
static unsigned int log2 (unsigned int val) {
    unsigned int ret = -1;
    while (val != 0) {
        val >>= 1;
        ret++;
    }
    return ret;
}
int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        

        // Extra options:
        //
        //IntOption    opt_verb   ("MAIN", "opt_verb",   "opt_verb level (0=silent, 1=some, 2=more).", 1, IntRange(0, 3));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
        


        StringOption    opt_graph("GRAPH", "graph","Not currently used", "");

        StringOption    opt_assume("MAIN", "assume","Specify a file of assumptions, with one literal or symbol per line", "");

        StringOption    opt_decidable("MAIN", "decidable-theories","Specify which graphs should make decisions on their own, in comma delimited format", "");

        StringOption 		opt_symbols("MAIN","symbols","Whether to read symbol lines (\"c var <variable number> <name>\") from the gnf","");

        StringOption    opt_assume_symbols("MAIN","assume-symbols","read in symbols (in the format produced by the 'symbols' option) and treat them as assumptions","");

        BoolOption opt_id_graph("GRAPH","print-vars","Identify the variables in the graph, then quit\n",false);

        BoolOption opt_witness("MAIN","witness","print solution",false);

        BoolOption   pre    ("MAIN", "pre",    "Completely turn on/off any preprocessing.", true);

        BoolOption opb("PB","opb","Parse the input as pseudo-boolean constraints in .opb format",false);

        parseOptions(argc, argv, true);

        bool using_symbols = strlen((const char* )opt_symbols)>0;

        if(opt_csv){
           	opt_verb=0;
           }
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        if(opt_verb>0)
        	fprintf(stderr,"WARNING: for repeatability, setting FPU to use double precision\n");
#endif

        vec<std::pair<int,string> > symbols;

        mincutalg = MinCutAlg::ALG_EDMONSKARP;

        if(!strcasecmp(opt_min_cut_alg,"ibfs")){
        	mincutalg=MinCutAlg::ALG_IBFS;

        }else if (!strcasecmp(opt_min_cut_alg,"edmondskarp-adj")){
        	mincutalg = MinCutAlg::ALG_EDKARP_ADJ;
        }else if (!strcasecmp(opt_min_cut_alg,"edmondskarp")){
        	mincutalg = MinCutAlg::ALG_EDMONSKARP;
        }else if  (!strcasecmp(opt_min_cut_alg,"edmondskarp-dynamic")){
        	mincutalg = MinCutAlg::ALG_EDKARP_DYN;
        }else if  (!strcasecmp(opt_min_cut_alg,"dinics")){
        	//Dinitz is also commonly spelled 'dinics' or 'Dinits', so accept those too...
        	mincutalg = MinCutAlg::ALG_DINITZ;
        }else if  (!strcasecmp(opt_min_cut_alg,"dinics-linkcut")){
        	mincutalg = MinCutAlg::ALG_DINITZ_LINKCUT;
        }else if  (!strcasecmp(opt_min_cut_alg,"dinitz")){
        	mincutalg = MinCutAlg::ALG_DINITZ;
        }else if  (!strcasecmp(opt_min_cut_alg,"dinitz-linkcut")){
        	mincutalg = MinCutAlg::ALG_DINITZ_LINKCUT;
        }else if  (!strcasecmp(opt_min_cut_alg,"dinits")){
        	//Dinitz is also commonly spelled 'dinics' or 'Dinits', so accept those too...
        	mincutalg = MinCutAlg::ALG_DINITZ;
        }else if  (!strcasecmp(opt_min_cut_alg,"dinits-linkcut")){
        	mincutalg = MinCutAlg::ALG_DINITZ_LINKCUT;
        }else{
        	fprintf(stderr,"Error: unknown max-flow/min-cut algorithm %s, aborting\n",((string)  opt_min_cut_alg).c_str());
        	exit(1);
        }

       componentsalg =ComponentsAlg::ALG_DISJOINT_SETS;

		 if(!strcasecmp(opt_components_alg,"disjoint-sets")){
			componentsalg=ComponentsAlg::ALG_DISJOINT_SETS;
		 }else{
			fprintf(stderr,"Error: unknown connectivity algorithm %s, aborting\n", ((string) opt_components_alg).c_str());
			exit(1);
		 }
		mstalg = MinSpanAlg::ALG_KRUSKAL;

		 if(!strcasecmp(opt_mst_alg,"kruskal")){
			 mstalg=MinSpanAlg::ALG_KRUSKAL;
		 }else if(!strcasecmp(opt_mst_alg,"prim")){
			 mstalg=MinSpanAlg::ALG_PRIM;
		 }else if (!strcasecmp(opt_mst_alg,"spira-pan")){
			 mstalg = MinSpanAlg::ALG_SPIRA_PAN;
		 }else{
			fprintf(stderr,"Error: unknown minimum spanning tree algorithm %s, aborting\n", ((string) opt_mst_alg).c_str());
			exit(1);
		 }

        reachalg = ReachAlg::ALG_BFS;

		 if(!strcasecmp(opt_reach_alg,"dijkstra")){
			reachalg=ReachAlg::ALG_DIJKSTRA;

		 }else if(!strcasecmp(opt_reach_alg,"bfs")){
			reachalg=ReachAlg::ALG_BFS;

		 }else if (!strcasecmp(opt_reach_alg,"dfs")){
			 reachalg = ReachAlg::ALG_DFS;
		 }else if (!strcasecmp(opt_reach_alg,"sat")){
			 reachalg = ReachAlg::ALG_SAT;
		 }else if (!strcasecmp(opt_reach_alg,"ramal-reps")){
			 reachalg = ReachAlg::ALG_RAMAL_REPS;
		 }else{
			fprintf(stderr,"Error: unknown reachability algorithm %s, aborting\n", ((string) opt_reach_alg).c_str());
			exit(1);
		 }


		 distalg = DistAlg::ALG_DISTANCE;

		 if(!strcasecmp(opt_dist_alg,"dijkstra")){
			 distalg=DistAlg::ALG_DIJKSTRA;

		 }else if(!strcasecmp(opt_dist_alg,"bfs")){
			 distalg=DistAlg::ALG_DISTANCE;

		 }else if (!strcasecmp(opt_dist_alg,"sat")){
			 distalg = DistAlg::ALG_SAT;
		 }else if (!strcasecmp(opt_dist_alg,"ramal-reps")){
			 distalg = DistAlg::ALG_RAMAL_REPS;
		 }else{
			fprintf(stderr,"Error: unknown distance algorithm %s, aborting\n", ((string) opt_dist_alg).c_str());
			exit(1);
		 }



		 undirectedalg = ConnectivityAlg::ALG_BFS;

		 if(!strcasecmp(opt_con_alg,"dijkstra")){
			 undirectedalg=ConnectivityAlg::ALG_DIJKSTRA;

		 }else if(!strcasecmp(opt_con_alg,"bfs")){
			 undirectedalg=ConnectivityAlg::ALG_BFS;

		 }else if (!strcasecmp(opt_con_alg,"dfs")){
			 undirectedalg = ConnectivityAlg::ALG_DFS;
		 }else if (!strcasecmp(opt_con_alg,"sat")){
			 undirectedalg = ConnectivityAlg::ALG_SAT;
		 }else  if (!strcasecmp(opt_con_alg,"thorup")){
			 undirectedalg = ConnectivityAlg::ALG_THORUP;
		 }/*else if (!strcasecmp(opt_con_alg,"ramal-reps")){
			 undirectedalg = ConnectivityAlg::ALG_RAMAL_REPS;
		 } */else{
			fprintf(stderr,"Error: unknown undirected reachability algorithm %s, aborting\n", ((string) opt_reach_alg).c_str());
			exit(1);
		 }

		    allpairsalg = AllPairsAlg::ALG_DIJKSTRA_ALLPAIRS;

		    if (!strcasecmp(opt_allpairs_alg,"floyd-warshall")){
		    		allpairsalg = AllPairsAlg::ALG_FLOYDWARSHALL;
		   		 }else if(!strcasecmp(opt_allpairs_alg,"dijkstra")){
		   			allpairsalg=AllPairsAlg::ALG_DIJKSTRA_ALLPAIRS;

			 } else{
				fprintf(stderr,"Error: unknown allpairs reachability algorithm %s, aborting\n", ((string) opt_allpairs_alg).c_str());
				exit(1);
			 }


		    undirected_allpairsalg = AllPairsConnectivityAlg::ALG_DIJKSTRA_ALLPAIRS;

			if (!strcasecmp(opt_undir_allpairs_alg,"floyd-warshall")){
				undirected_allpairsalg = AllPairsConnectivityAlg::ALG_FLOYDWARSHALL;
				 }else if(!strcasecmp(opt_undir_allpairs_alg,"dijkstra")){
					undirected_allpairsalg=AllPairsConnectivityAlg::ALG_DIJKSTRA_ALLPAIRS;
			 }else  if (!strcasecmp(opt_undir_allpairs_alg,"thorup")){
				 undirected_allpairsalg = AllPairsConnectivityAlg::ALG_THORUP;
			 } else{
				fprintf(stderr,"Error: unknown undirected allpairs reachability algorithm %s, aborting\n", ((string) opt_allpairs_alg).c_str());
				exit(1);
			 }

        double initial_time = cpuTime();

        const char * graphstr = opt_graph;


        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
#if not defined(__MINGW32__)
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    fprintf(stderr,"WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    fprintf(stderr,"WARNING! Could not set resource limit: Virtual memory.\n");
            } }
#endif
         const char *error;
         SimpSolver S;
         solver = &S;
         if (!pre) S.eliminate(true);
         S.max_decision_var = opt_restrict_decisions;
#ifdef DEBUG_SOLVER
         S.dbg_solver = new Solver();S.dbg_solver.verbosity=0;
#endif
         gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
             if (in == NULL)
                 printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

             if (opt_verb > 0){
                 printf("============================[ Problem Statistics ]=============================\n");
                 printf("|                                                                             |\n"); }

             Dimacs<StreamBuffer,SimpSolver> parser;
             GraphParser<char *, SimpSolver> graphParser;
             graphParser.setSymbols(&symbols);
             parser.addParser(&graphParser);
             parser.parse_DIMACS(in,S);

             gzclose(in);

             if(opt_verb>2){
            	 for(int i = 0;i<symbols.size();i++){
            		 int v = symbols[i].first;
            		 string s = symbols[i].second;
            		 std::cout<<"Symbol: " << (v+1) << " = " << s <<"\n";
            	 }
             }

         // Change to signal-handlers that will only notify the solver and allow it to terminate
           // voluntarily:
#if not defined(__MINGW32__)
           signal(SIGINT, SIGINT_interrupt);
           signal(SIGXCPU,SIGINT_interrupt);
#endif
     //      printf("Solving circuit with %d gates, %d latches, %d inputs, %d outputs\n", aiger->num_ands, aiger->num_latches, aiger->num_inputs, aiger->num_outputs);

           const char * priority_file = opt_priority;
           if(strlen(priority_file)>0){
        	   FILE * f = fopen(priority_file,"r");
        	   if(f){
        		   char * line = NULL;
				  int v=0;
				  int p=0;
				  int total_read = 0;
				  while (fscanf(f," %d %d ", &v,&p) ==2) {
					  if(v<1 || v> S.nVars() || p<0){
						  fprintf(stderr,"Bad priority line: %d %d", v, p);
						  exit(1);
					  }
					  v--;
					  total_read++;
					  S.setDecisionPriority(v,p);
					 }
				  if(total_read==0){
					  fprintf(stderr,"Warning: read no priorities from priority file!");
				  }
        		   fclose(f);
        	   }else{
        		   fprintf(stderr,"Failed to read priority file!\n");
        	   }
           }

           vec<int> decidable;

           string dstr= (const char*) opt_decidable;
           if(dstr.length()>0){
           std::replace(dstr.begin(),dstr.end(), '\'',' ');
           std::replace(dstr.begin(),dstr.end(), '\"',' ');
           std::replace(dstr.begin(),dstr.end(), ',',' ');

           istringstream iss(dstr);

            do
            {
                string sub;
                iss >> sub;
                if(sub.length()==0)
                	continue;
				//int value = atoi(sub.c_str());

                const char * s = sub.c_str();
                char* p=NULL;
                long value = strtol(s, &p, 10);
                if (*p) {
                    // conversion failed because the input wasn't a number
                }
                else {

                	 decidable.push(value);
                }

            } while (iss);
           }else{
        	   //default to all theories decidable
        	   for(int i = 0;i<S.theories.size();i++){
        		   decidable.push(i);
        	   }
           }
         //really simple, unsophisticated incremental BMC:
           vec<Lit> assume;

           const char * assume_str =opt_assume;
		   if(strlen(assume_str)){
	/*		 //gzFile gin =gzopen(assume_str, "rb");
			   int fd = open(assume_str, O_RDONLY, 0);
			   gzFile gin =   gzdopen(fd,"rb");
			parse_Assumptions(gin,assume, S,&symbols);
			if(assume.size()==0){
				printf("Warning: no assumptions found in assume file %s\n", assume_str);
			}
			gzclose(gin);*/
	     	   //FILE * f = fopen(assume_str,"r");
			   std::ifstream infile(assume_str);
			   std::string symbol;


			    std::unordered_map<std::string, int> symbol_table;
			    if(using_symbols){
			    	for (int i = 0;i<symbols.size();i++){
			    		int var = symbols[i].first;
			    		string symbol = symbols[i].second;
			    		symbol_table[symbol]=var;
			    	}
			    }

			   while (std::getline(infile, symbol))
			   {
				   std::stringstream trimmer;
				   trimmer << symbol;
				   symbol.clear();
				   trimmer>>symbol;
				   if(symbol.length()>0){
					   bool neg=false;
					   if(symbol[0]=='-'){
						   neg=true;
						   symbol.erase(0,1);
					   }
					   if(symbol.size()>0){
							if(symbol_table.count(symbol)==0){
								printf("PARSE ERROR! Unknown symbol: %s\n", symbol.c_str()), exit(3);
							}else{
								int v = symbol_table[symbol];
								Lit l =mkLit(v,neg);
								assume.push(l);
								if(opt_verb>2){
									if(neg)
										printf("Assume not %s (%d)\n", symbol.c_str(),dimacs(l));
									else
										printf("Assume %s = (%d)\n", symbol.c_str(),dimacs(l));
								}
							}
						}else{
							printf("PARSE ERROR! Empty symbol!\n"), exit(3);
						}
				   }
			   }
			    if(assume.size()==0){
			    	printf("Warning: no assumptions read from %s\n",assume_str);
			    }
		   }
		   if(opt_verb>2 && assume.size()){

			   printf("Assumptions: ");
				 for(int i = 0;i<assume.size();i++){
					 Lit l = assume[i];
					 printf("%d, ", dimacs(l));
				 }
				printf("\n");
		}

		   if(opt_verb>0 && decidable.size()){
			   printf("Decidable theories: ");
		   }
		   for(int i = 0;i<decidable.size();i++){
			   int t = decidable[i];
			   if(t<0 || t>= S.theories.size()){
				   fprintf(stderr,"Cannot set theory %d to be decidable, because there is no such theory\n",t);
				   fflush(stderr);
				   exit(1);
			   }
			   if(opt_verb>0){
					   printf("%d, ", t);
				   }
			   S.decidable_theories.push(S.theories[t]);
		   }
		   if(opt_verb>0 && decidable.size()){
				   printf("\n");
			   }

		   if(opt_id_graph){
			   if(S.theories.size()){

					   printf("Graph Variables:\n");
						Theory * t = S.theories[0];
						GraphTheorySolver *g = (GraphTheorySolver*)t;
						int width = sqrt(g->nNodes());
						if(opt_width>0){
							width=opt_width;
						}
						int height =width;
						if(opt_height>0){
							height = opt_height;
						}

						int v = 0;
						//for (int i = 0;i<w;i++){
						//	for(int j = 0;j<w;j++){
						int lasty= 0;
						for(int n = 0;n<height*width;n++){
							int x = n%width;
							int y = n/width;
							if(y > lasty)
								printf("\n");


							printf("%4d", n+1);

							Lit l = mkLit(n,false);
							if(assume.contains(l)){
								printf("+");
							}else{
								l=~l;
								if(assume.contains(l)){
									printf("-");
								}else{
									printf(" ");
								}
							}

							if (x<width-1){
								printf(",");
							}
							lasty=y;
						}
						printf("\n\n");




				}else{
					printf("No graph to identify\n");
				}

			   exit(0);
		   }

		if(strlen((const char* )opt_assume_symbols)>0){

			std::ifstream infile((const char* )opt_assume_symbols);

			std::string line;

			std::unordered_map<string,int> symbolmap;
			for(auto p:symbols){
				symbolmap[p.second]=p.first;
			}

			while (std::getline(infile, line))
			{
				if(line[0]=='%')
					continue;

			    std::istringstream iss(line);

			    string s1;
			    string s2;
			    string symbol;
			    bool sign=true;
			    iss>>s1 >> s2;
			    if(s1.compare(":-")){
			    	assert(false);
			    	fprintf(stderr,"Bad assumption: %s\n",line.c_str());
			    	exit(1);
			    }else if(!s2.compare("not")){
			    	//Then this is a _true_ assumption (yes, its intentionally backward...)
			    	iss>>symbol;
			    	sign=false;
			    }else{
			    	symbol=s2;
			    }
			    if(symbol.back()!='.'){
			      	assert(false);
					fprintf(stderr,"Bad assumption: %s\n",line.c_str());
					exit(1);
			    }
			    symbol.pop_back();
			    if(!symbolmap.count(symbol)){
			    	assert(false);
					fprintf(stderr,"Unmapped assumption symbol: %s\n",symbol.c_str());
					exit(1);
			    }
			    int var = symbolmap[symbol];
			    assert(var>=0);

			    Lit a = mkLit(var,sign);
			    assume.push(a);

			    /*line.
			    if(!line.compare(0, "2", ":-")){
			    	iss>>
			    }else{
			    	assert(false);
			    }*/


			}
		}

		double before_pre_processing = rtime(0);
		printf("simplify:\n");
		fflush(stdout);
		if(pre)
			S.eliminate(true);
		fflush(stdout);
		//exit(0);
		double preprocessing_time = rtime(0)-before_pre_processing;
		if(opt_verb>0 && pre){
			printf("Preprocessing time = %f\n", preprocessing_time);
		}
        lbool ret=S.solve(assume)?l_True:l_False;
        if(opt_optimize_mst && ret ==l_True && S.theories.size()){

        	GraphTheorySolver * g = (GraphTheorySolver*)S.theories[0];
        	if(g->mstDetector){

        		int mst_weight = g->mstDetector->positive_reach_detector->weight();

        		Var prev_min = var_Undef;
        		while(true){
        			printf("Optimizing minimum spanning tree (%d)...\n",mst_weight);
        			if(mst_weight==588983){
        				int a=1;
        			}
        			S.cancelUntil(0);
					Var min = S.newVar(true,false);
					if(min==3921){
						int a=1;
					}
					g->minimumSpanningTree(min,mst_weight-1);
					assume.push(mkLit(min,false));
					if(!S.solve(assume)){
						assume.pop();
						assume.push(mkLit(prev_min,false));
						bool check = S.solve(assume);
						assert(check);
						if(!check){
							fprintf(stderr,"Major Error! Instance is no longer satisfiable after removed assertition!\n");
							exit(3);
						}
						if(opt_check_solution){
											if(!g->check_solved()){
												fprintf(stderr,"Error! Solution doesn't satisfy theory properties!\n");
												exit(3);
											}
										}
						break;
					}
					mst_weight = g->mstDetector->positive_reach_detector->weight();
					assume.pop();
					prev_min = min;
					if(opt_check_solution){
										if(!g->check_solved()){
											fprintf(stderr,"Error! Solution doesn't satisfy theory properties!\n");
											exit(3);
										}
									}
        		}

        	}
		}

        if(ret==l_True){
        	if(!opt_csv)
        		printf("s SATISFIABLE\n");
        	 for(int i = 0;i<S.theories.size();i++)
        		 S.theories[i]->printSolution();

			if(opt_witness){

				printf("v ");
				for(int v =0;v<S.nVars();v++){
					if(S.model[v]==l_True){
						printf("%d ",(v+1));
					}else if(S.model[v]==l_False){
						printf("%d ",-(v+1));
					}/*else{
						printf("%d,",0);
					}*/
				}
				printf("0\n");
			}

			if(using_symbols){
				FILE * sfile = fopen(opt_symbols,"w");
				fprintf(sfile,"%% Generated by monosat\n");
				for(auto p:symbols){
					Var v = p.first;
					string & s = p.second;
					if(S.model[v]==l_True){
						fprintf(sfile,":- not %s.\n",s.c_str());
						//cout<<":- not "<< s<<".\n";
					}else if (S.model[v]==l_False){
						fprintf(sfile,":- %s.\n",s.c_str());
						//cout<<":- "<<s<<".\n";
					}else{
						//this is unassigned
						int a =1;
					}
				}
				fflush(sfile);
				fclose(sfile);
			}

			if(opb){
				printf("v ");
				for(auto p:symbols){
					Var v = p.first;
					std::string & s = p.second;
					if(S.model[v]==l_True){
						printf("%s ",s.c_str());
						//cout<<":- not "<< s<<".\n";
					}else if (S.model[v]==l_False){
						printf("-%s ",s.c_str());
						//cout<<":- "<<s<<".\n";
					}
				}
				printf("\n");
			}

        }else if(ret==l_False){
        	printf("s UNSATISFIABLE\n");
        }else{
        	printf("UNKNOWN\n");
        }
		if(opt_verb>0){
			printStats(S);

		}
        fflush(stdout);
#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}
