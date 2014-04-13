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
using namespace Minisat;

//=================================================================================================


void printStats(Solver& solver)
{
    double cpu_time = cpuTime();

    double mem_used = memUsedPeak();

    printf("restarts              : %" PRIu64 "\n", solver.starts);
    printf("conflicts             : %-12" PRIu64 "   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    printf("decisions             : %-12" PRIu64 "   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    printf("propagations          : %-12" PRIu64 "   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
    printf("conflict literals     : %-12" PRIu64 "   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
    printf("CPU time              : %g s\n", cpu_time);

    solver.printStats();
}


static Solver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt();   fflush(stdout);}

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
    if (opt_verb > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
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
        
        IntOption    opt_width("GRAPH","width","Width of graph.\n", 0, IntRange(0, INT32_MAX));
        IntOption    opt_height("GRAPH","height","Height of graph.\n", 0, IntRange(0, INT32_MAX));
        IntOption    opt_bits("GRAPH","bits","Bits per position in graph.\n", 1, IntRange(0, INT32_MAX));

        BoolOption	 opt_csv("GRAPH","csv","Output in CSV format",false);

        StringOption    opt_graph("GRAPH", "graph","Not currently used", "");

        StringOption    opt_assume("MAIN", "assume","Specify a file of assumptions, with one literal or symbol per line", "");

        StringOption    opt_decidable("MAIN", "decidable-theories","Specify which graphs should make decisions on their own, in comma delimited format", "");

        BoolOption 		opt_symbols("MAIN","read-symbols","Whether to read symbol lines (\"c var <variable number> <name>\") from the gnf",true);

        BoolOption opt_id_graph("GRAPH","print-vars","Identify the variables in the graph, then quit\n",false);

        BoolOption opt_witness("MAIN","witness","print solution",false);

        BoolOption   pre    ("MAIN", "pre",    "Completely turn on/off any preprocessing.", true);

        parseOptions(argc, argv, true);

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


        reachalg = ReachAlg::ALG_BFS;

		 if(!strcasecmp(opt_reach_alg,"dijkstra")){
			reachalg=ReachAlg::ALG_DIJKSTRA;

		 }else if(!strcasecmp(opt_reach_alg,"bfs")){
			reachalg=ReachAlg::ALG_BFS;

		 }else if (!strcasecmp(opt_reach_alg,"dfs")){
			 reachalg = ReachAlg::ALG_DFS;
		 }else if (!strcasecmp(opt_reach_alg,"sat")){
			 reachalg = ReachAlg::ALG_SAT;
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
		 }else{
			fprintf(stderr,"Error: unknown distance algorithm %s, aborting\n", ((string) opt_reach_alg).c_str());
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
		 } else{
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
         if (!pre) S.eliminate(true);
         S.max_decision_var = opt_restrict_decisions;
#ifdef DEBUG_SOLVER
         S.dbg_solver = new Solver();
#endif
         gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
             if (in == NULL)
                 printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

             if (opt_verb > 0){
                 printf("============================[ Problem Statistics ]=============================\n");
                 printf("|                                                                             |\n"); }

             parse_GRAPH(in, S,opt_symbols?&symbols:NULL);
             gzclose(in);

             if(opt_verb>2){
            	 for(int i = 0;i<symbols.size();i++){
            		 int v = symbols[i].first;
            		 string s = symbols[i].second;
            		 std::cout<<"Symbol: " << (v+1) << " = " << s <<"\n";
            	 }
             }

             if(strlen(graphstr)){
				 gzFile gin =gzopen(graphstr, "rb");
				parse_GRAPH(gin,S);
				gzclose(gin);
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
			    if(opt_symbols){
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

		   if(opt_verb>0){
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
		   if(opt_verb>0){
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
		S.eliminate(true);
        lbool ret=S.solve(assume)?l_True:l_False;
        if(opt_optimize_mst && ret ==l_True){

        	GraphTheorySolver * g = (GraphTheorySolver*)S.theories[0];
        	if(g->mstDetector){

        		int mst_weight = g->mstDetector->positive_reach_detector->weight();

        		Var prev_min = var_Undef;
        		while(true){
        			printf("Optimizing minimum spanning tree (%d)...\n",mst_weight);
        			if(mst_weight==6){
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
						if(opt_check_solution){
											if(!g->check_solved()){
												fprintf(stderr,"Error! Solution doesn't satisfy graph properties!\n");
												exit(1);
											}
										}
						break;
					}
					mst_weight = g->mstDetector->positive_reach_detector->weight();
					assume.pop();
					prev_min = min;
					if(opt_check_solution){
										if(!g->check_solved()){
											fprintf(stderr,"Error! Solution doesn't satisfy graph properties!\n");
											exit(1);
										}
									}
        		}

        	}
		}

        if(ret==l_True){
        	if(!opt_csv)
        		printf("s SATISFIABLE\n");

        	if(S.theories.size()){
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
				int bits = 1;
				if(opt_bits>0)
						bits=opt_bits;
				int v = 0;
				//for (int i = 0;i<w;i++){
				//	for(int j = 0;j<w;j++){
				int lasty= 0;
				int maxwidth = log10(pow(2, bits))+1; //highestbit(bits);
				for(int n = 0;n<height*width*bits;n+=bits){
					int x = n%(width*bits)/bits;
					int y = n/(width*bits);
					if(y > lasty)
						printf("\n");
#if not defined(__MINGW32__)
						if (!opt_csv && isatty(fileno(stdout))){
#else
						if(false){
#endif
							unsigned long val = 0;
							for(int j = 0;j<bits;j++){
								if(S.model[n+j]==l_True){
									val = val + (1<<j);
								}
							}

							//if(val>0){
								int backcolor = 0;
								if(val>0){
									backcolor=log2(val)+1;
								}
								if(backcolor<0){
									int a=1;
								}
								int forecolor = 7;
								if(backcolor>7){
									backcolor=7;
								}
								if(backcolor==3 || backcolor==7){
									forecolor=0;
								}
								printf("\033[1;4%dm\033[1;3%dm%*lu \033[0m",backcolor,forecolor,maxwidth,val);
							//}else{
								//printf("\033[1;44m\033[1;37m%*lu \033[0m",maxwidth,val);
								//printf("\033[1;40m\033[1;30m%*lu \033[0m",maxwidth,val);
							//}
						}else if (opt_csv){
							unsigned long val = 0;
							for(int j = 0;j<bits;j++){
								if(S.model[n+j]==l_True){
									val = val + (1<<j);
								}
							}
							printf("%*lu",maxwidth,val);
							if (x<width-1){
								printf(",");
							}
						}else{
							unsigned long val = 0;
							for(int j = 0;j<bits;j++){
								if(S.model[n+j]==l_True){
									val = val + (1<<j);
								}
							}
							printf(" %*lu ",maxwidth,val);

				/*			if(S.model[n]==l_True)
								printf(" 1");
							else
								printf(" 0");*/
						}

					lasty=y;
				}
				printf("\n\n");
				if(opt_check_solution){
							if(!g->check_solved()){
								fprintf(stderr,"Error! Solution doesn't satisfy graph properties!\n");
								exit(1);
							}
						}

				if(opt_print_reach){
				 v = 0;
				//for (int i = 0;i<w;i++){
				//	for(int j = 0;j<w;j++){
			/*	 lasty= 0;
				for(int n = 0;n<g->nNodes();n++){
					int x = n%width;
					int y = n/width;
					if(y > lasty)
						printf("\n");
#if not defined(__MINGW32__)
						if (isatty(fileno(stdout))){
#else
						if(false){
#endif

							if(S.model[n]==l_True)
								printf("\033[1;42m\033[1;37m%3d\033[0m",n);
							else
								printf("\033[1;44m\033[1;37m%3d\033[0m",n);
						}else{

							if(S.model[n]==l_True)
								printf(" 1");
							else
								printf(" 0");
						}

					lasty=y;
				}

				printf("\n");printf("\n");
				*/


				for(int t = 0;t<S.theories.size();t++){
					printf("Theory %d\n", t);
					GraphTheorySolver *g = (GraphTheorySolver*)S.theories[t];
					int nnodes = g->nNodes();

					int maxw = log10(g->nNodes() )+1; //highestbit(bits);

					{

						for(int r = 0;r<g->reach_detectors.size();r++){

							int width = sqrt(g->nNodes());
							if(opt_width>0){
									width=opt_width;
								}
								int height =width;
								if(opt_height>0){
									height = opt_height;
								}
							int lasty= 0;
							int extra =  g->nNodes() % width ? (width- g->nNodes() % width ):0;
							for(int n = 0;n<g->nNodes();n++){
								int x = n%width;

								int y = (n + extra )/width;
								if(y > lasty)
									printf("\n");

								int v =var( g->reach_detectors[r]->reach_lits[n]);
#if not defined(__MINGW32__)
								if (isatty(fileno(stdout)))
#else
								if(false)
#endif
								{
										if(S.model[v]==l_True)
											printf("\033[1;42m\033[1;37m%4d\033[0m", v+1);
										else
											printf("\033[1;44m\033[1;37m%4d\033[0m",v+1);
									}else{

										if(S.model[v]==l_True)
											printf(" 1");
										else
											printf(" 0");
									}

									lasty=y;
								}
								printf("\n");
							}



							//g->drawFull();

							assert(g->dbg_solved());
						}

					{
								for(int r = 0;r<g->distance_detectors.size();r++){

											int width = sqrt(g->nNodes());
											if(opt_width>0){
													width=opt_width;
												}
												int height =width;
												if(opt_height>0){
													height = opt_height;
												}
											int lasty= 0;
											int extra =  g->nNodes() % width ? (width- g->nNodes() % width ):0;
											for(int n = 0;n<g->nNodes();n++){
												int x = n%width;

												int y = (n + extra )/width;
												if(y > lasty)
													printf("\n");

												int d = g->distance_detectors[r]->positive_reach_detector->distance(n);
												printf("%*d ",maxw,d);


													lasty=y;
												}
												printf("\n");
											}
							}
						if(g->mstDetector){
								int min_weight = g->mstDetector->positive_reach_detector->weight();
								printf("Min Spanning Tree Weight: %d\n",min_weight);
								int width = sqrt(g->nNodes());
								if(opt_width>0){
										width=opt_width;
									}
									int height =width;
									if(opt_height>0){
										height = opt_height;
									}
								int lasty= 0;
								vec<bool> down_edge;
								int extra =  g->nNodes() % width ? (width- g->nNodes() % width ):0;
								for(int n = 0;n<g->nNodes();n++){
									int x = n%width;

									int y = (n + extra )/width;
									if(y > lasty){
										printf("\n");

										for(int i = 0;i<down_edge.size();i++){
											if(down_edge[i]){
												printf("|");
											}else{
												printf(" ");
											}
											printf(" ");
										}
										down_edge.clear();
										printf("\n");
									}
									printf("*");
									if(x<width-1){
										int edge_left = g->getEdgeID(n,n+1);
										Var edge_var = g->edge_list[edge_left].v;
										if(S.value(edge_var)==l_True &&  g->mstDetector->positive_reach_detector->edgeInTree(edge_left)){
											printf("-");
										}else{
											printf(" ");
										}
									}

									if(y<height-1){
											int edge_down = g->getEdgeID(n,n+width);
											Var edge_var = g->edge_list[edge_down].v;
											bool in_tree = g->mstDetector->positive_reach_detector->edgeInTree(edge_down);
											if(S.value(edge_var)==l_True &&  in_tree){
												down_edge.push(true);
											}else{
												down_edge.push(false);
											}
										}

										lasty=y;
									}
									printf("\n");
								}

						if(g->component_detector){
							int numComponents = g->component_detector->positive_reach_detector->numComponents();
							printf("Number of connected components is: %d\n",numComponents);

						}
					}




				}
/*        		for(int r = 0;r<g->reach_detectors.size();r++){

					int width = sqrt(g->nNodes());
					int lasty= 0;
					int extra =  g->nNodes() % width ? (width- g->nNodes() % width ):0;
					for(int n = 0;n<g->nNodes();n++){
						int x = n%width;

						int y = (n + extra )/width;

						int v =var( g->reach_detectors[r]->reach_lits[n]);
						if(v==306){
							int a =1;
						}
						if(S.model[v]==l_True){
							assert(S.value(v)==l_True);
							int node = g->reach_detectors[r]->getNode(v);
							g->reach_detectors[r]->positive_reach_detector->dbg_path(node);
							int  b=1;

						}


						lasty=y;
					}
					printf("\n");
				}*/
        	}

			if(opt_witness){

				printf("v ");
				for(int v =0;v<S.nVars();v++){
					if(S.model[v]==l_True){
						printf("+%d ",(v+1));
					}else if(S.model[v]==l_False){
						printf("%d ",-(v+1));
					}/*else{
						printf("%d,",0);
					}*/
				}
				printf("0\n");
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
