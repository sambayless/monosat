/***************************************************************************************[Solver.cc]
 The MIT License (MIT)

 Copyright (c) 2018, Sam Bayless

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

package monosat;

import java.util.Optional;

/**
 * Literals are integers in the rance 0..nVars()*2
 * Literals come in pairs; positive literals are even, negative literals are odd.
 * Finally, there are two reserved literals:
 * lit_Undef and lit_Error, with values -2 and -1, respectively.
 * <p>
 * In an ideal world, this would be _just_ an integer, or a type-checked value class of size 32 bits.
 * However, to deal with
 */
public final class Lit {
    public final static Lit Undef = new Lit(-1, true);
    public final static Lit Error = new Lit(-2, true);
    public final static Lit True = new Lit(0, true);
    public final static Lit False = new Lit(1, true);
    protected final Solver solver;

    //The value of this literal
    protected int l = -2;

    private Lit(int lit, boolean define_literal) {
        //used to define static lits
        assert(lit<=1);
        this.solver=null;
        this.l = lit;
    }

    protected Lit(Solver solver,int literal) {
        assert (literal >= 0);
        this.l =literal;
        this.solver=solver;
    }

    /**
     * Before you can create a new Literal, you must first create a Solver:
     * Solver s = new Solver();
     * Lit a = new Lit(s);
     * @param solver The SAT solver in which to create a new literal.
     */
    public Lit(Solver solver){
        this(solver,true);
    }

    /**
     * Before you can create a new Literal, you must first create a Solver:
     * Solver s = new Solver();
     * Lit a = new Lit(s);
     * @param solver The SAT solver in which to create a new literal.
     * @param decisionLit If true (the default), this literal can be used as a decision Literal in the solver.
     */
    public Lit(Solver solver, boolean decisionLit){
        int var = MonosatJNI.newVar(solver.solverPtr);
        this.l = var*2;
        this.solver = solver;
        solver.registerLit(this);
        solver.setDecisionLiteral(this,decisionLit);
    }

    public Solver getSolver(){
        if(solver==null){
            throw new RuntimeException("Cannot access solver for " + toString());
        }
        return solver;
    }

    /**
     * Check whether the literal has a negation sign.
     *
     * @return True if the literal is negative, false if it is positive.
     */
    public boolean sign() {
        return (l & 1) == 1;
    }

    public Lit not() {
        if(this==True){
            return False;
        }else if (this==False){
            return True;
        }
        return solver.not(this);
    }
    public Lit abs() {
        if (sign()) {
            return this.not();
        } else {
            return this;
        }
    }

    /**
     * Convert the literal into dimacs format
     *
     * @return
     */
    public int toDimacs() {
        return ((l / 2) + 1) * (this.sign() ? -1 : 1);
    }

    @Override
    public String toString() {
        if(l==-2){
            return "Lit_Error";
        }else if(l==-1){
            return "Lit_Undef";
        }else if (l==0){
            return "True";
        }else if (l==1){
            return "False";
        }
        return "Lit" + l;
    }

    /**
     * Convert the literal into a variable.
     *
     * @return
     */
    public int toVar() {
        return l / 2;
    }


    /**
     * Convert the literal into an integer.
     *
     * @return
     */
    public int toInt() {
        return l;
    }

    /**
     * Query the model in the solver, throwing a NoModelException if the literal is unassigned in the model.
     * This can happen only if the literal is not a decision literal.
     */
    public boolean value() throws NoModelException {
        return value(Solver.LBool.Undef);
    }
    /**
     * Query the model in the solver.
     * If defaultVal is LBool.Undef, this will throw a NoModelException if the literal is unassigned.
     * Else, if the literal is unassigned, defaultVal will be returned.
     */
    public boolean value(boolean defaultVal) throws NoModelException {
        return value(Solver.LBool.fromBool(defaultVal));
    }
    /**
     * After a solve call, non-decision literals may or may not be assigned to a value.
     * Unassigned literals will have the value LBool.Undef;
     */
    public Optional<Boolean> possibleValue(){
        return getPossibleValue(Solver.LBool.Undef);
    }
    public Optional<Boolean> possibleValue(boolean defaultValue){
        return getPossibleValue(defaultValue ? Solver.LBool.True : Solver.LBool.False);
    }

    /**
     * After a solve call, non-decision literals may or may not be assigned to a value.
     * Unassigned literals will have the value defaultValue
     */
    private Optional<Boolean> getPossibleValue( Solver.LBool defaultValue) {
        if(this.l<0){
            throw new IllegalArgumentException("Literal " + toString() + " is not a valid literal.");
        }
        if(this==Lit.True){
            return Solver.LBool.True.toOpt();
        }else if(this==Lit.False){
            return Solver.LBool.False.toOpt();
        }

        Solver.LBool val = Solver.LBool.toLbool(MonosatJNI.getModel_Literal(getSolver().solverPtr, l));
        if (val==Solver.LBool.Undef){
            return defaultValue.toOpt();
        }else{
            return val.toOpt();
        }
    }

    protected Solver.LBool getLBoolValue(Solver.LBool defaultValue) {

        Solver.LBool val = Solver.LBool.toLbool(MonosatJNI.getModel_Literal(getSolver().solverPtr, l));
        if (val==Solver.LBool.Undef){
            return defaultValue;
        }else{
            return val;
        }
    }


    /**
     * Query the model in the solver.
     * If defaultVal is LBool.Undef, this will throw an exception if the literal is unassigned.
     * Else, if the literal is unassigned, defaultVal will be returned.
     */
    protected boolean value(Solver.LBool defaultVal) throws NoModelException {
        if(this.l<0){
            throw new IllegalArgumentException("Literal " + toString() + " is not a valid literal.");
        }
        if(this==Lit.True){
            return true;
        }else if(this==Lit.False){
            return false;
        }

        Solver.LBool val = Solver.LBool.toLbool(MonosatJNI.getModel_Literal(getSolver().solverPtr, l));
        if (val == Solver.LBool.Undef) {
            if (defaultVal == Solver.LBool.Undef) {
                throw new NoModelException("Literal " + toString() + " is unassigned in the current model");
            } else {
                val = defaultVal;
            }
        }
        return val == Solver.LBool.True;
    }

    /**
     * Sometimes the solver may prove that a literal must always be true, or must always be false, in any satisfying assignment.
     * If the solver has already proven that literal l must always be true, or always false, returns LBool.True or LBool.False.
     * Otherwise, returns LBool.Undef.
     *
     * Note that as constraints are added to the solver, and after solve calls, the solver may learn that new literals are constant.
     * This function will not attempt to prove that l is constant, if the solver does not already know this.
     * You can attempt to discover whether a literal is constant using the assumption mechanism, eg:
     * if(~solve(l)){
     *   //l must be constant false, because the constraints are UNSAT if l is true.
     * }
     *
     */
    public Optional<Boolean> getConstantValue() {
        if(this.l<0){
            throw new IllegalArgumentException("Literal " + toString() + " is not a valid literal.");
        }
        if(this==Lit.True){
            return Solver.LBool.True.toOpt();
        }else if(this==Lit.False){
            return Solver.LBool.False.toOpt();
        }
        Solver.LBool val = Solver.LBool.toLbool(MonosatJNI.getConstantModel_Literal(solver.solverPtr, l));
        return val.toOpt();
    }


    /**
     * If a literal is known to the solver to be a constant (either true or false),
     * this returns True.
     * Note that even if this function returns false, the literal may still in fact be a constant,
     * but this may not yet be known to the solver. As constraints are added, or after calls to solve(),
     * the solver may discover that a literal is constant, and so this return value may change over time.
     * @return True if the literal is known to the solver to be either always true, or always false, in any satisfying model.
     */
    public boolean isConst(){
        if(this==Lit.True || this==Lit.False){
            return true;
        }else if(l<0){
            return false;
        }
        return getConstantValue().isPresent();
    }

    public boolean isConstTrue(){
        if(this==Lit.True){
            return true;
        }else if (this==Lit.False || this==Lit.Undef || this==Lit.Error){
            return false;
        }
        return getConstantValue().orElse(false);
    }
    public boolean isConstFalse(){
        if(this==Lit.False){
            return true;
        }else if (this==Lit.True || this==Lit.Undef || this==Lit.Error){
            return false;
        }
        return !getConstantValue().orElse(true);
    }
}
