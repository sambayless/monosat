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
import monosat.Logic;
import monosat.Solver;
/**
 * Literals are integers in the rance 0..nVars()*2
 * Literals come in pairs; positive literals are even, negative literals are odd.
 * Finally, there are two reserved literals:
 * lit_Undef and lit_Error, with values -2 and -1, respectively.
 * <p>
 * In an ideal world, this would be _just_ an integer, or a type-checked value class of size 32 bits.
 * However, to deal with
 */
final public class Lit {
    public final static Lit Undef = new Lit(-1, true);
    public final static Lit Error = new Lit(-2, true);
    public final static Lit True = new Lit(0, true);
    public final static Lit False = new Lit(1, true);
    protected final Solver solver;

    //private final Lit neg; //every literal also has a pointer to its negation.
    //The value of this literal
    protected int l = -2;

/*
    *//**
     * Typically, users will create literals using Solver.newLit(), rather than constructing them directly.
     *//*
    protected Lit(Solver solver) {
        this.l = -2;
        this.solver=solver;
        //this.neg = Error;
    }*/

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
     * Query the model in the solver, throwing an exception if the literal is unassigned in the model.
     * This can happen only if the literal is not a decision literal.
     */
    public boolean value() throws RuntimeException {
        if(this==Lit.True) {
            return true;
        }else if (this==Lit.False){
            return false;
        }
        return getSolver().getValue(this);
    }
    /**
     * Query the model in the solver.
     * If defaultVal is LBool.Undef, this will throw an exception if the literal is unassigned.
     * Else, if the literal is unassigned, defaultVal will be returned.
     */
    public boolean value(LBool defaultVal) throws RuntimeException {
        if(this==Lit.True) {
            return true;
        }else if (this==Lit.False){
            return false;
        }
        return getSolver().getValue(this, defaultVal);
    }
    /**
     * After a solve call, non-decision literals may or may not be assigned to a value.
     * Unassigned literals will have the value LBool.Undef;
     */
    public LBool possibleValue(){
        if(this==Lit.True) {
            return LBool.True;
        }else if (this==Lit.False){
            return LBool.False;
        }
        return getSolver().getPossibleValue(this);
    }
    /**
     * If a literal is known to the solver to be a constant (either true or false),
     * this returns that value. Otherwise, returns LBool.Undef
     */
    public LBool constantValue(){
        if(this==Lit.True) {
            return LBool.True;
        }else if (this==Lit.False){
            return LBool.False;
        }
        return getSolver().getConstantValue(this);
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
        return getSolver().getConstantValue(this)!=LBool.Undef;
    }

    public boolean isConstTrue(){
        if(this==Lit.True){
            return true;
        }else if (this==Lit.False || this==Lit.Undef || this==Lit.Error){
            return false;
        }
        return getSolver().getConstantValue(this)==LBool.True;
    }
    public boolean isConstFalse(){
        if(this==Lit.False){
            return true;
        }else if (this==Lit.True || this==Lit.Undef || this==Lit.Error){
            return false;
        }
        return getSolver().getConstantValue(this)==LBool.False;
    }
}
