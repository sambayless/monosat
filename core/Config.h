/*
 * Config.h
 *
 *  Created on: 2012-10-20
 *      Author: sam
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include "utils/Options.h"
namespace Minisat {
extern  DoubleOption  opt_var_decay;
extern  DoubleOption  opt_clause_decay ;
extern  DoubleOption  opt_random_var_freq ;
extern  DoubleOption  opt_random_seed;
extern  IntOption     opt_ccmin_mode  ;
extern  IntOption     opt_phase_saving ;
extern  BoolOption    opt_rnd_init_act ;
extern  BoolOption    opt_luby_restart ;
extern  IntOption     opt_restart_first ;
extern  DoubleOption  opt_restart_inc   ;
extern  DoubleOption  opt_garbage_frac  ;

extern BoolOption opt_partial;

extern  BoolOption opt_interpolate;
extern IntOption opt_eager_prop;
extern IntOption opt_subsearch;


extern  IntOption opt_allsat_vars;
extern BoolOption opt_allsat_first;
extern BoolOption opt_allsat;
}

#endif /* CONFIG_H_ */
