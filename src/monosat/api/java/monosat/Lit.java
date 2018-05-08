/*
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
 */

package monosat;

import java.util.Optional;

/**
 * Represents a signed Boolean literal in MonoSAT.
 * Internally, literals are integers in the rance 0..(nVars()-1)*2.
 * Literals come in pairs; positive literals are even, negative literals are odd.
 * <br>
 * There are four constant literals that are globally shared among all Solver instances:
 * Lit.True and Lit.False, with values 0 and 1, and represent constant true and false, and <br>
 * Lit.Undef and Lit.Error, with values -1 and -2, respectively.
 * <br>
 * All other literals belong to specific Solver instances.<p>
 * Example usage:
 *
 * <blockquote><pre>{@code
 * Solver solver = new Solver();
 * Lit a = new Lit(solver);//create a new, free literal
 * Lit b = new Lit(solver);//create another literal
 *
 * solver.addClause(a,b.not()); //add a clause requiring either a to be true, or b to be false.
 * solver.addClause(a.not,b); //add a clause requiring either a to be false, or b to be true.
 * //Together, the above clauses for a to equal b in any satisfying solution.
 * //See also monosat.Logic.assertEqual(a,b).
 * solver.solve();//attempt to solve this formula
 *
 * //retrieve the values of each literal from the solver's satisfying assignment
 * bool a_assign = a.value();
 * bool b_assign = b.value();
 * }</pre></blockquote>
 */
public final class Lit {
    /**
     * This constant value represents an undefined literal.
     */
    public final static Lit Undef = new Lit(-1, true);
    /**
     * This constant value represents an invalid literal.
     */
    public final static Lit Error = new Lit(-2, true);
    /**
     * A constant literal, shared among all solver instances, that is true in any satisfying assignment.
     */
    public final static Lit True = new Lit(0, true);

    /**
     * A constant literal, shared among all solver instances, that is false in any satisfying assignment.
     */
    public final static Lit False = new Lit(1, true);

    /**
     * The solver instance this literal belongs to.
     */
    protected final Solver solver;

    /**
     * The integer value of this literal.
     * Note that this is not marked finall, as released literals may be renumbered.
     */
    protected int l = -2;

    /**
     * A private constructor used only to initialize the constant, global literals.
     * @param lit The integer value of this literal.
     * @param define_literal Ignored.
     */
    private Lit(int lit, boolean define_literal) {
        //used to define static lits
        assert(lit<=1);
        this.solver=null;
        this.l = lit;
    }

    /**
     * Used internally by the Solver to create literals.
     * @param solver Solver this literal will belong to.
     * @param literal The integer value of this literal.
     */
    protected Lit(Solver solver,int literal) {
        assert (literal >= 0);
        this.l =literal;
        this.solver=solver;
    }

    /**
     * Create a new literal in the specified solver.
     * Note: Before you can create a new Literal, you must first create a Solver:
     *
     * <blockquote><pre>{@code
     * Solver s = new Solver();
     * Lit a = new Lit(s);
     * }</pre></blockquote>
     *
     * @param solver The SAT solver in which to create a new literal.
     */
    public Lit(Solver solver){
        this(solver,true);
    }

    /**
     * Create a new literal in the specified solver.
     * Note: Before you can create a new Literal, you must first create a Solver:
     *
     * <blockquote><pre>{@code
     * Solver s = new Solver();
     * Lit a = new Lit(s);
     * }</pre></blockquote>
     *
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

    /**
     * Get the solver this literal belongs to.
     * Note that the 4 constant, global literals (Lit.True,Lit.False,Lit.Undef,Lit.Error)
     * do not have a solver.
     * @return The solver this literal belongs to.
     */
    public Solver getSolver(){
        if(solver==null){
            throw new RuntimeException("Cannot access solver for " + toString());
        }
        return solver;
    }

    /**
     * Check whether this literal has a negation sign.
     *
     * @return True if the literal is negative, false if it is positive.
     */
    public boolean sign() {
        return (l & 1) == 1;
    }

    /**
     * Return the negation of this literal.
     * <p>
     * Implementation note:<br>
     * Lit.not() is well optimized and inexpensive, in both Java and the native library.
     * Double negations will be eliminated and are essentially free.
     * No new objects are instantiated by Lit.not().
     * It is always the case that <code>a.not().not() == a<code/>
     *
     * @return The negation of this literal.
     */
    public Lit not() {
        if(this==True){
            return False;
        }else if (this==False){
            return True;
        }
        if(this==Error || this==Undef){
            throw new RuntimeException("Cannot negate literal " + toString());
        }
        return solver.not(this);
    }

    /**
     * Return this.not() if this has a negation sign, else return this.
     * @return this.not() if this has a negation sign, else return this.
     */
    public Lit abs() {
        if (sign()) {
            return this.not();
        } else {
            return this;
        }
    }

    /**
     * Convert the literal into DIMACS format.
     * Dimacs format literals are non-zero, and are negative if the literal has negative sign.
     * For example, Lit.True.toDimacs() = 1, and Lit.False.toDimacs() = -1.
     *
     * @return An integer representation of the literal in dimacs format.
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
     * @return An integer representation of the variable this literal represents (ignoring sign).
     */
    public int toVar() {
        return l / 2;
    }


    /**
     * Convert the literal into an integer.
     *
     * @return An integer representation of this literal, even if the literal is positive, odd if it is negative.
     */
    public int toInt() {
        return l;
    }

    /**
     * Query the model in the solver, throwing a NoModelException if the literal is unassigned in the model.
     * This can happen only if the literal is not a decision literal.
     *
     * @return The value of this literal in the solvers model (if the solver has a satisfying assignment.)
     * @throws NoModelException If the solver does not yet have a satisfying assignment (eg, if s.solve() has not yet
     * returned true.)
     */
    public boolean value() throws NoModelException {
        return value(Solver.LBool.Undef);
    }

    /**
     * Query the model in the solver.     *
     * Else, if the literal is unassigned, defaultVal will be returned.
     *
     * @param defaultVal If the literal is unassigned in the solver, returns defaultVal. Else, returns the value of
     *                   the literal in the solver.
     * @return The value of this literal in the solvers model (if the solver has an assignment for this
     * literal), else defaultVal.
     * @throws NoModelException If the solver does not yet have a satisfying assignment (eg, if s.solve() has not yet
     * returned true.)
     */
    public boolean value(boolean defaultVal) {
        return value(Solver.LBool.fromBool(defaultVal));
    }

    /**
     * After a solve call, non-decision literals may or may not be assigned to a value.
     * Unassigned literals will have the value LBool.Undef;
     *
     * @return An Optional containing the value of this literal if it has an assignment.
     */
    public Optional<Boolean> possibleValue(){
        return getPossibleValue(Solver.LBool.Undef);
    }

    /**
     * After a solve call, non-decision literals may or may not be assigned to a value.
     * Unassigned literals will have the value defaultValue
     *
     * @return An Optional containing the value of this literal if it has an assignment, else return  defaultValue
     * .toOpt().
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

    /**
     * Query the model in the solver.
     * If defaultVal is LBool.Undef, this will throw an exception if the literal is unassigned.
     * Else, if the literal is unassigned, defaultVal will be returned.
     *
     * @param defaultVal The value to return if this literal is undefined.
     * @return The value of this literal in the model if it has a value, else if defaultValue != LBool.Undef,
     * return defaultValue==Solver.LBool.True; else, throw a NoModelException
     * @throws NoModelException If there is no assignment for this literal in the solver,
     * and defaultVal is LBool.Undef.
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
     * @return An optional containing the constant value of this literal if the solver knows it to be trivially true
     * or false; Optional.empty if the solver has not proven that the literal is  constant.
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
     *
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


    /**
     * If a literal is known to the solver to be a constant (either true or false),
     * this returns True.
     * Note that even if this function returns false, the literal may still in fact be a constant,
     * but this may not yet be known to the solver. As constraints are added, or after calls to solve(),
     * the solver may discover that a literal is constant, and so this return value may change over time.
     *
     * @return True if the literal is known to the solver to be either always true in any satisfying model.
     */
    public boolean isConstTrue(){
        if(this==Lit.True){
            return true;
        }else if (this==Lit.False || this==Lit.Undef || this==Lit.Error){
            return false;
        }
        return getConstantValue().orElse(false);
    }


    /**
     * If a literal is known to the solver to be a constant (either true or false),
     * this returns True.
     * Note that even if this function returns false, the literal may still in fact be a constant,
     * but this may not yet be known to the solver. As constraints are added, or after calls to solve(),
     * the solver may discover that a literal is constant, and so this return value may change over time.
     *
     * @return True if the literal is known to the solver to be either always false in any satisfying model.
     */
    public boolean isConstFalse(){
        if(this==Lit.False){
            return true;
        }else if (this==Lit.True || this==Lit.Undef || this==Lit.Error){
            return false;
        }
        return !getConstantValue().orElse(true);
    }
}
