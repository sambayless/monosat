//
// Created by sam on 06/07/16.
//

#ifndef MONOSAT_CONFIG_PB_H
#define MONOSAT_CONFIG_PB_H
#include "minisatpb/ADTs/Int.h"
#include "minisatpb/ADTs/Global.h"

namespace Monosat {
namespace PB {
enum SolverT {
    st_MiniSat, st_SatELite
};
enum ConvertT {
    ct_Sorters, ct_Adders, ct_BDDs, ct_Mixed, ct_Undef
};
enum Command {
    cmd_Minimize, cmd_FirstSolution, cmd_AllSolutions
};

// -- output options:
//these two moved to Global.h
/*extern bool     opt_satlive;
extern bool     opt_ansi;*/
extern char *opt_cnf;
extern int opt_verbosity;
extern bool opt_try;

// -- solver options:
extern ConvertT opt_convert;
extern ConvertT opt_convert_goal;
extern bool opt_convert_weak;
extern double opt_bdd_thres;
extern double opt_sort_thres;
extern double opt_goal_bias;
extern Int opt_goal;
extern Command opt_command;
extern bool opt_branch_pbvars;
extern int opt_polarity_sug;

// -- files:
extern char *opt_input;
extern char *opt_result;
}
}
#endif //MONOSAT_CONFIG_PB_H
