//
// Created by sam on 06/07/16.
//

#ifndef MONOSAT_CONFIG_PB_H
#define MONOSAT_CONFIG_PB_H

#include "monosat/pb/ADTs/Int.h"
#include "monosat/pb/ADTs/Global.h"
#include "monosat/utils/Options.h"


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

extern IntOption opt_verbosity;
extern BoolOption opt_try;

// -- solver options:
extern ConvertT opt_convert;
extern ConvertT opt_convert_goal;
extern BoolOption opt_convert_weak;
extern DoubleOption opt_bdd_thres;
extern DoubleOption opt_sort_thres;
extern DoubleOption opt_goal_bias;
extern Int64Option opt_goal;
extern Command opt_command;
extern BoolOption opt_branch_pbvars;
extern IntOption opt_polarity_sug;
extern BoolOption opt_preprocess;
}
}
#endif //MONOSAT_CONFIG_PB_H
