//
// Created by sam on 12/07/16.
//

#ifndef MONOSAT_PB_H
#define MONOSAT_PB_H
#include "monosat/core/SolverTypes.h"
namespace Monosat {
namespace PB {
enum Ineq{
    LT = -2,
    LEQ = -1,
    EQ = 0,
    GEQ = 1,
    GT = 2
};
//abstract interface to Psuedoboolean constraint solver
class PBConstraintSolver {
public:
    virtual bool addConstr(const vec<Lit> &ps, const vec<int> &Cs, int rhs, Ineq ineq)=0;
    virtual Lit addConditionalConstr(const vec<Lit> &ps, const vec<int> &Cs, int rhs, Ineq ineq, Lit cond = lit_Undef)=0;

};
}
}
#endif //MONOSAT_PB_H
