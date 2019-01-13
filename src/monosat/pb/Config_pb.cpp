#include "Config_pb.h"

using namespace Monosat::PB;
static const char* _cat_pb = "PB";

Monosat::BoolOption  Monosat::PB::opt_try(_cat_pb, "pb-try",
                                          "(hidden option -- if set, then \"try\" to parse, but don't output \"s UNKNOWN\" if you fail, instead exit with error code 5)",
                                          false);

Monosat::IntOption Monosat::PB::opt_verbosity(_cat_pb, "pb-verb",
                                              "", 0, IntRange(0, 3));
Monosat::BoolOption Monosat::PB::opt_preprocess(_cat_pb, "pb-pre", "", true);
Monosat::BoolOption Monosat::PB::opt_convert_weak(_cat_pb, "pb-convert-weak", "", true);
Monosat::DoubleOption Monosat::PB::opt_bdd_thres(_cat_pb, "pb-bdd-threshold", "", 3,
                                                 DoubleRange(0, false, HUGE_VAL, false));
Monosat::DoubleOption Monosat::PB::opt_sort_thres(_cat_pb, "pb-sort-threshold", "", 20,
                                                  DoubleRange(0, false, HUGE_VAL, false));
Monosat::DoubleOption Monosat::PB::opt_goal_bias(_cat_pb, "pb-goal-bias", "", 3,
                                                 DoubleRange(0, false, HUGE_VAL, false));
Monosat::Int64Option Monosat::PB::opt_goal(_cat_pb, "pb-goal-bias", "", 9223372036854775807LL,
                                           Int64Range(INT64_MIN, INT64_MAX));
Monosat::BoolOption Monosat::PB::opt_branch_pbvars(_cat_pb, "pb-branch-pbvars", "", false);
Monosat::IntOption Monosat::PB::opt_polarity_sug(_cat_pb, "pb-polarity", "", 1, IntRange(-1, 1));

ConvertT Monosat::PB::opt_convert = ct_Mixed;
ConvertT Monosat::PB::opt_convert_goal = ct_Undef;
Command Monosat::PB::opt_command = cmd_Minimize;

