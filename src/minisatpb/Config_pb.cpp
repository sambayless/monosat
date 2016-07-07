#include "minisatpb/Config_pb.h"


char*    opt_cnf       = NULL;
int      opt_verbosity = 1;
bool     opt_try       = false;     // (hidden option -- if set, then "try" to parse, but don't output "s UNKNOWN" if you fail, instead exit with error code 5)

bool     opt_preprocess    = true;
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
                "MiniSat+ 1.1, based on MiniSat 2.2.0  -- (C) Niklas Een, Niklas Sorensson, 2012\n"
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
                "USAGE: minisatp <input-file> [<result-file>] [-<option> ...]\n"
                "\n"
                "Solver options:\n"
                "  -ca -adders   Convert PB-constrs to clauses through adders.\n"
                "  -cs -sorters  Convert PB-constrs to clauses through sorters.\n"
                "  -cb -bdds     Convert PB-constrs to clauses through bdds.\n"
                "  -cm -mixed    Convert PB-constrs to clauses by a mix of the above. (default)\n"
                "  -ga/gs/gb/gm  Override conversion for goal function (long name: -goal-xxx).\n"
                "  -w -weak-off  Clausify with equivalences instead of implications.\n"
                "  -no-pre       Don't use MiniSat's CNF-level preprocessing.\n"
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

            else if (oneof(arg, "ca,adders" )) opt_convert = ct_Adders;
            else if (oneof(arg, "cs,sorters")) opt_convert = ct_Sorters;
            else if (oneof(arg, "cb,bdds"   )) opt_convert = ct_BDDs;
            else if (oneof(arg, "cm,mixed"  )) opt_convert = ct_Mixed;

            else if (oneof(arg, "ga,goal-adders" )) opt_convert_goal = ct_Adders;
            else if (oneof(arg, "gs,goal-sorters")) opt_convert_goal = ct_Sorters;
            else if (oneof(arg, "gb,goal-bdds"   )) opt_convert_goal = ct_BDDs;
            else if (oneof(arg, "gm,goal-mixed"  )) opt_convert_goal = ct_Mixed;

            else if (oneof(arg, "w,weak-off"     )) opt_convert_weak = false;
            else if (oneof(arg, "no-pre"))          opt_preprocess   = false;

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