/*****************************************************************************[PbSolver_convert.cc]
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

#include "PbSolver.h"
#include "Hardware.h"

namespace Monosat {
namespace PB {
//-------------------------------------------------------------------------------------------------
void linearAddition(const Linear& c, vec<Formula>& out);        // From: PbSolver_convertAdd.C
Formula buildConstraint(const Linear& c, int max_cost = INT_MAX);   // From: PbSolver_convertSort.C
Formula convertToBdd(const Linear& c, int max_cost = INT_MAX);   // From: PbSolver_convertBdd.C
//-------------------------------------------------------------------------------------------------


bool PbSolver::convertPbs(bool first_call){
    vec<Formula> converted_constrs;

    if(first_call){
        findIntervals();
        if(!rewriteAlmostClauses()){
            sat_solver.addEmptyClause();
            return false;
        }
    }

    for(int i = 0; i < constrs.size(); i++){
        if(constrs[i] == NULL) continue;
        Linear& c = *constrs[i];
        assert(c.lo != Int_MIN || c.hi != Int_MAX);

        if(opt_verbosity >= 1)
            /**/reportf("---[%4d]---> ", constrs.size() - 1 - i);

        if(opt_convert == ct_Sorters)
            converted_constrs.push(buildConstraint(c));
        else if(opt_convert == ct_Adders)
            linearAddition(c, converted_constrs);
        else if(opt_convert == ct_BDDs)
            converted_constrs.push(convertToBdd(c));
        else if(opt_convert == ct_Mixed){
            int adder_cost = estimatedAdderCost(c);
            //**/printf("estimatedAdderCost: %d\n", estimatedAdderCost(c));
            Formula result = convertToBdd(c, (int) (adder_cost * opt_bdd_thres));
            if(result == _undef_)
                result = buildConstraint(c, (int) (adder_cost * opt_sort_thres));
            if(result == _undef_)
                linearAddition(c, converted_constrs);
            else
                converted_constrs.push(result);
        }else
            assert(false);

        if(!okay()) return false;
    }

    // NOTE: probably still leaks memory (if there are constraints that are NULL'ed elsewhere)
    for(int i = 0; i < constrs.size(); i++)
        if(constrs[i] != NULL)
            constrs[i]->~Linear();

    constrs.clear();
    mem.clear();

    clausify(*this, converted_constrs);

    return okay();
}
}
}