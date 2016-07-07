/*****************************************************************************************[Main.cc]
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


#include <unistd.h>
#include <signal.h>
#include "PbSolver.h"
#include "PbParser.h"
#include "Config_pb.h"

using namespace Monosat::PB;

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


static void SIGINT_handler(int /*signum*/) {
    reportf("\n");
    reportf("*** INTERRUPTED ***\n");
    //SatELite::deleteTmpFiles();
    _exit(0); }     // (using 'exit()' rather than '_exit()' sometimes causes the solver to hang (why?))


static void SIGTERM_handler(int /*signum*/) {
    reportf("\n");
    reportf("*** TERMINATED ***\n");
    outputResult(*pb_solver, false);
    //SatELite::deleteTmpFiles();
    _exit(0);
}


PbSolver::solve_Command convert(Command cmd) {
    switch (cmd){
    case cmd_Minimize:      return PbSolver::sc_Minimize;
    case cmd_FirstSolution: return PbSolver::sc_FirstSolution;
    default: 
        assert(cmd == cmd_AllSolutions);
        return PbSolver::sc_AllSolutions;
    }
}


//=================================================================================================


int main(int argc, char** argv)
{
    /*DEBUG*/if (argc > 1 && (strcmp(argv[1], "-debug") == 0 || strcmp(argv[1], "--debug") == 0)){ void test(); test(); exit(0); }

    parseOptions(argc, argv);
    pb_solver = new PbSolver(opt_preprocess);
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
        pb_solver->printStats();
        reportf("_______________________________________________________________________________\n");
    }

    exit(0); // (faster than "return", which will invoke the destructor for 'PbSolver')
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

    SimpSolver      S;
    vec<Formula>    fs;
    fs.push(f ^ g);
    clausify(S, fs);

    S.verbosity = 1;
    printf(S.solve() ? "SAT\n" : "UNSAT\n");
}
