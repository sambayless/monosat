/******************************************************************************************[Main.C]
Copyright (c) 2005-2010, Niklas Een, Niklas Sorensson

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

/**************************************************************************************************

Read a DIMACS file and apply the SAT-solver to it.

**************************************************************************************************/


#include <cstdarg>
#include <unistd.h>
#include <signal.h>
#include "MiniSat.h"
#include "PbSolver.h"
#include "PbParser.h"
#include "Main.h"

//=================================================================================================
// Command line options:


bool     opt_satlive   = true;
bool     opt_ansi      = true;
char*    opt_cnf       = NULL;
int      opt_verbosity = 1;
bool     opt_try       = false;     // (hidden option -- if set, then "try" to parse, but don't output "s UNKNOWN" if you fail, instead exit with error code 5)

SolverT  opt_solver        = st_MiniSat;
ConvertT opt_convert       = ct_Mixed;
ConvertT opt_convert_goal  = ct_Undef;
bool     opt_convert_weak  = true;
double   opt_bdd_thres     = 3;
double   opt_sort_thres    = 20;
double   opt_goal_bias     = 3;
Int      opt_goal          = Int_MAX;
Command  opt_command       = cmd_Minimize;
bool     opt_branch_pbvars = false;
int      opt_polarity_sug  = 1;
bool     opt_old_format    = false;

char*    opt_input  = NULL;
char*    opt_result = NULL;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

cchar* doc =
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    "MiniSat+ 1.0, based on MiniSat v1.13  -- (C) Niklas Een, Niklas Sorensson, 2005\n"
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    "USAGE: minisat+ <input-file> [<result-file>] [-<option> ...]\n"
    "\n"
    "Solver options:\n"
    "  -M -minisat   Use MiniSat v1.13 as backend (default)\n"
    "  -S -satelite  Use SatELite v1.0 as backend\n"
    "\n"
    "  -ca -adders   Convert PB-constrs to clauses through adders.\n"
    "  -cs -sorters  Convert PB-constrs to clauses through sorters.\n"
    "  -cb -bdds     Convert PB-constrs to clauses through bdds.\n"
    "  -cm -mixed    Convert PB-constrs to clauses by a mix of the above. (default)\n"
    "  -ga/gs/gb/gm  Override conversion for goal function (long name: -goal-xxx).\n"
    "  -w -weak-off  Clausify with equivalences instead of implications.\n"
    "\n"
    "  -bdd-thres=   Threshold for prefering BDDs in mixed mode.        [def: %g]\n"
    "  -sort-thres=  Threshold for prefering sorters. Tried after BDDs. [def: %g]\n"
    "  -goal-bias=   Bias goal function convertion towards sorters.     [def: %g]\n"
    "\n"
    "  -1 -first     Don\'t minimize, just give first solution found\n"
    "  -A -all       Don\'t minimize, give all solutions\n"
    "  -goal=<num>   Set initial goal limit to '<= num'.\n"
    "\n"
    "  -p -pbvars    Restrict decision heuristic of SAT to original PB variables.\n"
    "  -ps{+,-,0}    Polarity suggestion in SAT towards/away from goal (or neutral).\n"
    "\n"
    "Input options:\n"
    "  -of -old-fmt  Use old variant of OPB file format.\n"
    "\n"
    "Output options:\n"
    "  -s -satlive   Turn off SAT competition output.\n"
    "  -a -ansi      Turn off ANSI codes in output.\n"
    "  -v0,-v1,-v2   Set verbosity level (1 default)\n"
    "  -cnf=<file>   Write SAT problem to a file. Trivial UNSAT => no file written.\n"
    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

bool oneof(cchar* arg, cchar* alternatives)
{
    // Force one leading '-', allow for two:
    if (*arg != '-') return false;
    arg++;
    if (*arg == '-') arg++;

    // Scan alternatives:
    vec<char*>  alts;
    splitString(alternatives, ",", alts);
    for (int i = 0; i < alts.size(); i++){
        if (strcmp(arg, alts[i]) == 0){
            xfreeAll(alts);
            return true;
        }
    }
    xfreeAll(alts);
    return false;
}


void parseOptions(int argc, char** argv)
{
    vec<char*>  args;   // Non-options

    for (int i = 1; i < argc; i++){
        char*   arg = argv[i];
        if (arg[0] == '-'){
            if (oneof(arg,"h,help")) fprintf(stderr, doc, opt_bdd_thres, opt_sort_thres, opt_goal_bias), exit(0);

            else if (oneof(arg, "M,minisat" )) opt_solver = st_MiniSat;
            else if (oneof(arg, "S,satelite")) opt_solver = st_SatELite;

            else if (oneof(arg, "ca,adders" )) opt_convert = ct_Adders;
            else if (oneof(arg, "cs,sorters")) opt_convert = ct_Sorters;
            else if (oneof(arg, "cb,bdds"   )) opt_convert = ct_BDDs;
            else if (oneof(arg, "cm,mixed"  )) opt_convert = ct_Mixed;

            else if (oneof(arg, "ga,goal-adders" )) opt_convert_goal = ct_Adders;
            else if (oneof(arg, "gs,goal-sorters")) opt_convert_goal = ct_Sorters;
            else if (oneof(arg, "gb,goal-bdds"   )) opt_convert_goal = ct_BDDs;
            else if (oneof(arg, "gm,goal-mixed"  )) opt_convert_goal = ct_Mixed;

            else if (oneof(arg, "w,weak-off"     )) opt_convert_weak = false;

            //(make nicer later)
            else if (strncmp(arg, "-bdd-thres=" , 11) == 0) opt_bdd_thres  = atof(arg+11);
            else if (strncmp(arg, "-sort-thres=", 12) == 0) opt_sort_thres = atof(arg+12);
            else if (strncmp(arg, "-goal-bias=",  11) == 0) opt_goal_bias  = atof(arg+11);
            else if (strncmp(arg, "-goal="     ,   6) == 0) opt_goal       = atoi(arg+ 6);  // <<== real bignum parsing here
            else if (strncmp(arg, "-cnf="      ,   5) == 0) opt_cnf        = arg + 5;
            //(end)

            else if (oneof(arg, "1,first"   )) opt_command = cmd_FirstSolution;
            else if (oneof(arg, "A,all"     )) opt_command = cmd_AllSolutions;

            else if (oneof(arg, "p,pbvars"  )) opt_branch_pbvars = true;
            else if (oneof(arg, "ps+"       )) opt_polarity_sug = +1;
            else if (oneof(arg, "ps-"       )) opt_polarity_sug = -1;
            else if (oneof(arg, "ps0"       )) opt_polarity_sug =  0;

            else if (oneof(arg, "of,old-fmt" )) opt_old_format = true;

            else if (oneof(arg, "s,satlive" )) opt_satlive = false;
            else if (oneof(arg, "a,ansi"    )) opt_ansi    = false;
            else if (oneof(arg, "try"       )) opt_try     = true;
            else if (oneof(arg, "v0"        )) opt_verbosity = 0;
            else if (oneof(arg, "v1"        )) opt_verbosity = 1;
            else if (oneof(arg, "v2"        )) opt_verbosity = 2;

            else
                fprintf(stderr, "ERROR! Invalid command line option: %s\n", argv[i]), exit(1);

        }else
            args.push(arg);
    }

    if (args.size() == 0)
        fprintf(stderr, doc, opt_bdd_thres, opt_sort_thres, opt_goal_bias), exit(0);
    if (args.size() >= 1)
        opt_input = args[0];
    if (args.size() == 2)
        opt_result = args[1];
    else if (args.size() > 2)
        fprintf(stderr, "ERROR! Too many files specified on commandline.\n"),
        exit(1);
}


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


void reportf(const char* format, ...)
{
    static bool col0 = true;
    static bool bold = false;
    va_list args;
    va_start(args, format);
    char* text = vnsprintf(format, args);
    va_end(args);

    for(char* p = text; *p != 0; p++){
        if (col0 && opt_satlive)
            putchar('c'), putchar(' ');

        if (*p == '\b'){
            bold = !bold;
            if (opt_ansi)
                putchar(27), putchar('['), putchar(bold?'1':'0'), putchar('m');
            col0 = false;
        }else{
            putchar(*p);
            col0 = (*p == '\n' || *p == '\r');
        }
    }
    fflush(stdout);
}


//=================================================================================================
// Helpers:


PbSolver*   pb_solver = NULL;   // Made global so that the SIGTERM handler can output best solution found.


void outputResult(const PbSolver& S, bool optimum = true)
{
    if (!opt_satlive) return;

    if (optimum){
        if (S.best_goalvalue == Int_MAX) printf("s UNSATISFIABLE\n");
        else                             printf("s OPTIMUM FOUND\n");
    }else{
        if (S.best_goalvalue == Int_MAX) printf("s UNKNOWN\n");
        else                             printf("s SATISFIABLE\n");
    }

    if (S.best_goalvalue != Int_MAX){
        printf("v");
        for (int i = 0; i < S.best_model.size(); i++)
            printf(" %s%s", S.best_model[i]?"":"-", S.index2name[i]);
        printf("\n");
    }
    fflush(stdout);
}


static void SIGINT_handler(int signum) {
    reportf("\n");
    reportf("*** INTERRUPTED ***\n");
    SatELite::deleteTmpFiles();
    _exit(0); }     // (using 'exit()' rather than '_exit()' sometimes causes the solver to hang (why?))


static void SIGTERM_handler(int signum) {
    reportf("\n");
    reportf("*** TERMINATED ***\n");
    outputResult(*pb_solver, false);
    SatELite::deleteTmpFiles();
    _exit(pb_solver->best_goalvalue == Int_MAX ? 0 : 10); }


void printStats(BasicSolverStats& stats, double cpu_time)
{
    reportf("restarts              : %"I64_fmt"\n", stats.starts);
    reportf("conflicts             : %-12"I64_fmt"   (%.0f /sec)\n", stats.conflicts   , stats.conflicts   /cpu_time);
    reportf("decisions             : %-12"I64_fmt"   (%.0f /sec)\n", stats.decisions   , stats.decisions   /cpu_time);
    reportf("propagations          : %-12"I64_fmt"   (%.0f /sec)\n", stats.propagations, stats.propagations/cpu_time);
    reportf("inspects              : %-12"I64_fmt"   (%.0f /sec)\n", stats.inspects    , stats.inspects    /cpu_time);
    reportf("CPU time              : %g s\n", cpu_time);
}


PbSolver::solve_Command convert(Command cmd) {
    switch (cmd){
    case cmd_Minimize:      return PbSolver::sc_Minimize;
    case cmd_FirstSolution: return PbSolver::sc_FirstSolution;
    case cmd_AllSolutions:  return PbSolver::sc_AllSolutions;
    default: assert(false); }
}


//=================================================================================================


int main(int argc, char** argv)
{
    /*DEBUG*/if (argc > 1 && (strcmp(argv[1], "-debug") == 0 || strcmp(argv[1], "--debug") == 0)){ void test(); test(); exit(0); }

    parseOptions(argc, argv);
    pb_solver = new PbSolver(); // (must be constructed AFTER parsing commandline options -- constructor uses 'opt_solver' to determinte which SAT solver to use)
    signal(SIGINT , SIGINT_handler);
    signal(SIGTERM, SIGTERM_handler);

    // Set command from 'PBSATISFIABILITYONLY':
    char* value = getenv("PBSATISFIABILITYONLY");
    if (value != NULL && atoi(value) == 1)
        reportf("Setting switch '-first' from environment variable 'PBSATISFIABILITYONLY'.\n"),
        opt_command = cmd_FirstSolution;

    if (opt_verbosity >= 1) reportf("Parsing PB file...\n");
    parse_PB_file(opt_input, *pb_solver, opt_old_format);


    pb_solver->solve(convert(opt_command));

    if (pb_solver->goal == NULL && pb_solver->best_goalvalue != Int_MAX)
        opt_command = cmd_FirstSolution;    // (otherwise output will be wrong)
    if (!pb_solver->okay())
        opt_command = cmd_Minimize;         // (HACK: Get "UNSATISFIABLE" as output)

    // <<== write result to file 'opt_result'

    if (opt_command == cmd_Minimize)
        outputResult(*pb_solver);
    else if (opt_command == cmd_FirstSolution)
        outputResult(*pb_solver, false);

    if (opt_verbosity >= 1) {
        reportf("_______________________________________________________________________________\n\n");
        printStats(pb_solver->stats, cpuTime());
        reportf("_______________________________________________________________________________\n");
    }

    exit(pb_solver->best_goalvalue == Int_MAX ? 20 : (pb_solver->goal == NULL || opt_command == cmd_FirstSolution) ? 10 : 30);    // (faster than "return", which will invoke the destructor for 'PbSolver')
}



//=================================================================================================
#include "Hardware.h"
#include "Debug.h"

#define N 10

void test(void)
{
    Formula f = var(0), g = var(N-1);
    for (int i = 1; i < N; i++)
        f = Bin_new(op_Equiv, f, var(i)),
        g = Bin_new(op_Equiv, g, var(N-i-1));

    dump(f); dump(g);

    printf("f= %d\n", index(f));
    printf("g= %d\n", index(g));

    Solver          S(true);
    vec<Formula>    fs;
    fs.push(f ^ g);
    clausify(S, fs);

    S.setVerbosity(1);
    printf(S.solve() ? "SAT\n" : "UNSAT\n");
}
