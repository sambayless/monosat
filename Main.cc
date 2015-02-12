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
#include <gmpxx.h>
#include <errno.h>
#include <climits>
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
#include "Cover.h"
#include "Subsume.h"
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

long n_blocking_clauses=0;
long total_clause_length=0;
mpz_class n_solutions=0;
//double n_solutions=0;
long max_clauses=0;
long max_int =0;
static Solver* solver;
FILE * intout;
void printInterpolant (Solver & S, FILE * out){
	fprintf(out,"c Over-approximation of satisfying assignments\n");
	for(int i = 0;i<S.interpolant.size();i++){
		vec<Lit> & c = S.interpolant[i];
		for(int j  =0;j<c.size();j++){
			Lit l = c[j];
			fprintf(out,"%d ",dimacs(l));
		}
		fprintf(out,"0\n");
	}
	fflush(out);
}


// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {

    if (solver->verbosity > 0){
        printStats(*solver);
        }
    printf("\n"); printf("*** INTERRUPTED ***\n");

    printf("Interrupted allsat\n");

		printf("Learned %ld blocking clauses.\n",n_blocking_clauses);
		printf("Avg. blocking clause length: %f\n",total_clause_length/((double)n_blocking_clauses));


   		if(opt_subsume){
   			Subsume subsume;
   			subsume.setNumVars(solver->nVars());
   			printf("Checking for subsumed clauses...\n");
   			subsume.checkAll(solver->interpolant);
   		}

   		long total_int_clause_length = 0;
   		int n_int_clauses =0;
   		for(int i = 0;i<solver->interpolant.size();i++){
   			n_int_clauses++;
   			total_int_clause_length+= solver->interpolant[i].size();
   		}

   		if(opt_interpolate){
   			printf("# Interpolants: %d\n", solver->interpolant.size());
   			printf("Avg. interpolant clause length: %f\n",total_int_clause_length/((double)n_int_clauses));
   		}
   		gmp_printf ("#solutions: %Zd\n", n_solutions.get_mpz_t());
	   if(intout){
			printInterpolant(*solver,intout);
			fclose(intout);
		}

	   fflush(stdout);


    _exit(1);
}


// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt(); fflush(stdout); SIGINT_exit(signum); }
//=================================================================================================
// Main:


int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        intout=NULL;
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    solver_verb   ("MAIN", "solver-verb",   "Verbosity level (0=silent, 1=some, 2=more).", 0, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
        
        StringOption opt_interpolant_file("MAIN","int-out","Write interpolant to given file","");

        parseOptions(argc, argv, true);
        if(!opt_allsat_modsat){
        	opt_interpolate=false;//interpolants are currently only supported by modsat
        }

        double initial_time = cpuTime();

        max_clauses= opt_max_allsat;
        max_int=opt_max_interpolant;
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


	   S.verbosity = solver_verb;
	   S.use_model=false;
	   solver = &S;
        if (S.verbosity > 0){
            printf("============================[ Problem Statistics ]=============================\n");
            printf("|                                                                             |\n"); }

        parse_DIMACS(in, S);
        gzclose(in);
        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;
        const char * intfile =  opt_interpolant_file;
        intout = strlen(intfile)>0? fopen(intfile,"wb") : NULL;

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
        bool has_any_solutions=false;
        if(!opt_allsat){

        	ret = S.solveLimited(dummy);
      	  if(ret==l_True) {
      		has_any_solutions=true;
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

        	for(int i = opt_allsat_from-1;i<S.nVars() && i < opt_allsat_to;i++){

				allsatvec.push(i);
			}
        	printf("Performing allsat on the following %d variables: ", allsatvec.size());
        	for (int i = 0;i<allsatvec.size();i++){
        		Var v = allsatvec[i];
        		printf("%d ", (v+1));
        	}
        	printf("\n");

        	//do an allsat loop
        	{
        		vec<Lit> learnt_clause;
        		vec<Var> supervec;
        		vec<Lit> block;
        		Cover c;
        		Solver allsat_solver;
        		Solver & allsat = opt_allsat_modsat?allsat_solver:S;
        		allsat.use_model=false;
        		for(int i = 0;i<S.nVars();i++){
					c.excludeFromCover(i,true);
				}
        		vec<int> allsat_map;


        		if(opt_allsat_modsat){
					//map from sub to supersolver vars

					for(int i = 0;i<allsatvec.size();i++){
						Var v =allsat.newVar();
						supervec.push(v);
						allsat_map.growTo(allsatvec[i]+1);
						allsat_map[allsatvec[i]]=v;
						c.excludeFromCover(allsatvec[i],false);
					}
					allsat.addTheory(&S);
					S.attachTo(&allsat,supervec,allsatvec);
				}else{
					for(int i = 0;i<allsatvec.size();i++){
						c.excludeFromCover(allsatvec[i],false);
						allsat_map.growTo(allsatvec[i]+1);
						allsat_map[allsatvec[i]]=allsatvec[i];
						if(opt_decide_allsat_first){
							allsat.setDecisionPriority(allsatvec[i],1);
						}
					}
				}

        		//ok, now do an allsat loop:
        		while(ret!=l_False && (!max_clauses || n_blocking_clauses<max_clauses) && (!max_int || S.interpolant.size() <max_int)  &&! allsat.asynch_interrupt){
        			ret = allsat.solveLimited(dummy);
        			if(ret==l_True){
        				has_any_solutions=true;
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
        					mpz_class p_solutions=0;
        					int width = (allsatvec.size() - block.size());
        					mpz_ui_pow_ui (p_solutions.get_mpz_t(), 2, width);
        					//double t = pow(2, );
        					n_solutions+= p_solutions;
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
        				if(opt_allsat_inc_block && allsat.decisionLevel()>0 && block.size()>0){
        					int max_lev = allsat.level(var(block[0]));
        					int second_max = block.size()>1 ? allsat.level(var(block[1])):0;
							int max_pos=0;
							int second_pos=1;
							Lit max=block[0];
							for(int i = 0;i<block.size();i++){
								Lit l = block[i];
								assert(allsat.value(l)==l_False);
								int lev = allsat.level(var(l));
								if(lev>=max_lev){
									second_max=max_lev;
									second_pos=max_pos;
									max_lev=lev;
									max_pos=i;
									max=l;
								}else if (lev>second_max){
									second_max=lev;
									second_pos=i;
								}
							}

							if(block.size()>1){
								Lit second_max = block[second_pos];
								if(max_pos!=0){
									Lit t = block[0];
									block[0]=max;
									block[max_pos]=t;
									if(second_pos==0){
										second_pos=max_pos;
									}
								}

								if(second_pos!=1){
									Lit t= block[1];
									block[1]=second_max;
									block[second_pos]=t;
								}

								assert(block[0]==max);
								assert(block[1]==second_max);
							}


							allsat.cancelUntil(second_max);

							if(block.size()==1 || max_lev>second_max){
								if (block.size() == 1){
									allsat.ok&= allsat.enqueue(max);
								}else{
									CRef cr = allsat.ca.alloc(block, false);
									allsat.clauses.push(cr);//previously, was adding these to blocking_clauses, which causes Cover to ignore them - which might possibly lead to double counting (not sure about this).
									allsat.attachClause(cr);

									assert(allsat.value(block[0])==l_Undef);
									allsat.uncheckedEnqueue(block[0], cr);
								}
							}else{
								//solver is in conflict...
								CRef confl = allsat.ca.alloc(block, false);
								allsat.clauses.push(confl);//blocking_clauses
								allsat.attachClause(confl);

								if(allsat.decisionLevel()==0){
									ret = l_False;//done
								}else{


									 learnt_clause.clear();
									 int backtrack_level=0;
									allsat.analyze(confl, learnt_clause, backtrack_level);
									allsat.cancelUntil(backtrack_level);

									//this is now slightly more complicated, if there are multiple lits implied by the super solver in the current decision level:
									//The learnt clause may not be asserting.

									if (learnt_clause.size() == 1){
										allsat.uncheckedEnqueue(learnt_clause[0]);
									}else{
										CRef cr = allsat.ca.alloc(learnt_clause, true);
										allsat.learnts.push(cr);
										allsat.attachClause(cr);
										allsat.claBumpActivity(allsat.ca[cr]);


										allsat.uncheckedEnqueue(learnt_clause[0], cr);



									}
								}
							}
        				}else{
        					allsat.cancelUntil(0);
        					allsat.addClause_(block);
        				}

        				total_clause_length+=block.size();

        			}
        		}

        		printf("Done allsat\n");

        		printf("Learned %ld blocking clauses.\n",n_blocking_clauses);
        		printf("Avg. blocking clause length: %f\n",total_clause_length/((double)n_blocking_clauses));
        		printf("# Decisions %d, # conflicts %d\n", allsat.decisions, allsat.conflicts);
        		if(opt_allsat_modsat)
        			printf("Subsolver: # Decisions %d, # conflicts %d\n", S.decisions, S.conflicts);
        		if(opt_interpolate){
					if(opt_subsume){
						Subsume subsume;
						subsume.setNumVars(S.nVars());
						printf("Checking for subsumed clauses...\n");
						subsume.checkAll(S.interpolant);
					}

					long total_int_clause_length = 0;
					int n_int_clauses =0;
					for(int i = 0;i<S.interpolant.size();i++){
						n_int_clauses++;
						total_int_clause_length+= S.interpolant[i].size();
					}

        			printf("# Interpolants: %d\n", S.interpolant.size());
        			printf("Avg. interpolant clause length: %f\n",total_int_clause_length/((double)n_int_clauses));
        		}
        		gmp_printf ("#solutions: %Zd\n", n_solutions.get_mpz_t());

       		   if(intout){
       		        printInterpolant(S,intout);
 				}
        	}
        }
        if(intout){
        	fclose(intout);
        }
        if(has_any_solutions){
        	ret=l_True;
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
