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

#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "Aiger.h"

#include "core/Solver.h"
#include "Aiger.h"
#include "core/Dimacs.h"
#include "core/Config.h"
#include "Partial.h"
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
    if (solver->verbosity > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
    fflush(stdout);
    _exit(1); }


//=================================================================================================
// Main:


int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
        

        parseOptions(argc, argv, true);


        double initial_time = cpuTime();



        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    printf("WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    printf("WARNING! Could not set resource limit: Virtual memory.\n");
            } }


        gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
        if (in == NULL)
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

        Solver S;


	   S.verbosity = verb;
	   S.use_model=false;
	   solver = &S;
        if (S.verbosity > 0){
            printf("============================[ Problem Statistics ]=============================\n");
            printf("|                                                                             |\n"); }

        parse_DIMACS(in, S);
        gzclose(in);
        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;

        if (S.verbosity > 0){
            printf("|  Number of variables:  %12d                                         |\n", S.nVars());
            printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }

        double parsed_time = cpuTime();
        if (S.verbosity > 0){
            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
            printf("|                                                                             |\n"); }

        // Change to signal-handlers that will only notify the solver and allow it to terminate
        // voluntarily:
        signal(SIGINT, SIGINT_interrupt);
        signal(SIGXCPU,SIGINT_interrupt);

        if (!S.simplify()){
            if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
            if (S.verbosity > 0){
                printf("===============================================================================\n");
                printf("Solved by unit propagation\n");
                printStats(S);
                printf("\n"); }
            printf("UNSATISFIABLE\n");
            exit(20);
        }
        vec<Lit> dummy;
        lbool ret=l_Undef;
        if(!opt_allsat){

        	ret = S.solveLimited(dummy);
      	  if(ret==l_True) {
      		 if(opt_partial){
      			Cover c;
      			vec<Lit> partial;
      			c.getCover(S,partial);
      			sort(partial);
 				printf("v ");
 				for (int i = 0; i < partial.size(); i++){
 					Lit l = partial[i];
 					printf("%s%s%d", (i==0)?"":" ", (sign(l))?"-":"", var(l)+1);

 				}
 				printf(" 0\n");


      		 }else{
				printf("v ");
				for (int i = 0; i < S.nVars(); i++){
					if (S.value(i) != l_Undef)
						printf("%s%s%d", (i==0)?"":" ", (S.value(i)==l_True)?"":"-", i+1);
				}

				printf(" 0\n");
      		 }
      	  }

        }else{
        	vec<Var> allsatvec;
        	if(opt_allsat_vars==0){
        		for(int i = 0;i<S.nVars();i++){

					allsatvec.push(i);
				}
        	}else if(opt_allsat_vars>0){
        		for(int i = 0;i<opt_allsat_vars;i++){
        			if(i>=S.nVars())
        				break;
        			allsatvec.push(i);
        		}
        	}else {
        		int start = S.nVars()+opt_allsat_vars;
        		if(start<0)
        			start=0;
        		for(int i = start;i<S.nVars();i++){
					allsatvec.push(i);
				}
        	}

        	//do an allsat loop
        	if(opt_allsat_first){
        		vec<Var> supervec;
        		vec<Lit> block;
        		Cover c;
        		Solver allsat;
        		allsat.use_model=false;
        		for(int i = 0;i<S.nVars();i++){
					c.excludeFromCover(i,true);
				}
        		//map from sub to supersolver vars
        		vec<int> allsat_map;

        		for(int i = 0;i<allsatvec.size();i++){
        			Var v =allsat.newVar();
        			supervec.push(v);
        			allsat_map.growTo(allsatvec[i]+1);
        			allsat_map[allsatvec[i]]=v;
        			c.excludeFromCover(allsatvec[i],false);
        		}
        		allsat.addTheory(&S);
        		S.attachTo(&allsat,supervec,allsatvec);
        		long n_blocking_clauses=0;
        		double n_solutions=0;
        		//ok, now do an allsat loop:
        		while(ret!=l_False){
        			ret = allsat.solveLimited(dummy);
        			if(ret==l_True){
        				n_blocking_clauses++;
        				//Negate the solver's assignments
        				block.clear();
        				if(opt_partial){
        					c.getCover(S,block);
        					sort(block);
        					if(verb>1){
								printf("v ");
								for (int i = 0; i < block.size(); i++){
									Lit l = block[i];
									printf("%s%s%d", (i==0)?"":" ", (sign(l))?"-":"", var(l)+1);

								}
								printf(" 0\n");
        					}
        					double t = pow(2, (allsatvec.size() - block.size()));
        					n_solutions+= t;
        					for(int i = 0;i<block.size();i++){
        						Lit l = block[i];
        						Var s =allsat_map[var(l)];
        						block[i]= ~mkLit(s,sign(l));//intentionally inverting this
        					}
        				}else{
        					n_solutions++;
        					for(int i = 0;i<allsat.nVars();i++){
								block.push(~mkLit(i, allsat.value(i)==l_False));
							}
          					if(verb>1){
								printf("v ");
								for (int i = 0; i < block.size(); i++){
									Lit l = block[i];
									l=mkLit(allsatvec[var(l)],sign(l));
									printf("%s%s%d", (i==0)?"":" ", (sign(l))?"-":"", var(l)+1);

								}
								printf(" 0\n");
							}

        				}
        				allsat.cancelUntil(0);
        				allsat.addClause(block);
        			}
        		}
        		printf("Done allsat\n");
        		printf("Learned %ld blocking clauses.\n",n_blocking_clauses);
        		printf("#solutions: %.0f\n",n_solutions);
        	}
        }
        if(ret==l_True){
        	printf("SAT\n");
        }else if(ret==l_False){
        	printf("UNSAT\n");
        }else{
        	printf("UNKNOWN\n");
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
