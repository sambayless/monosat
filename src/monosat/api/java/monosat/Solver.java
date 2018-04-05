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

import com.sun.xml.internal.bind.v2.model.util.ArrayInfoUtil;

import java.io.Closeable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

public class Solver implements Closeable {
    //Holds positive versions of all literals, so that we don't need to create multiple literal objects for the same literal
    //As literals have no reference to the solver and are really just a thin wrapper around an integer,
    //this list can be safely shared across all solvers.
    private static ArrayList<Lit> lits = new ArrayList<>();
    private static long MAX_CACHE_CONST = 255;
    protected long solverPtr = 0; //handle to the underlying monsoat this instance.
    protected long bvPtr = 0;
    int buffer_size0 = 1024;
    int buffer_size1 = 1024;
    int buffer_size2 = 1024;
    Lit true_lit;
    private IntBuffer ints0;
    private IntBuffer ints1;
    private IntBuffer ints2;
    private ArrayList<ArrayList<BitVector>> cached_bvs = new ArrayList<ArrayList<BitVector>>();

    /**
     * Instantiate a new Solver.
     * By default, support for preprocessing is disabled (as solving with preprocessing enabled requires some extra care)
     */
    public Solver() {
        this(false);
    }

    public Solver(String args) {
        this(args, false);
    }

    public Solver(boolean enable_preprocessing) {
        this("", enable_preprocessing);
    }

    public Solver(ArrayList<String> args) {
        this(collectArgs(args), false);
    }

    public Solver(ArrayList<String> args, boolean enable_preprocessing) {
        this(collectArgs(args), enable_preprocessing);
    }

    public Solver(String args, boolean enable_preprocessing) {
        if (args != null && args.length() > 0) {
            solverPtr = MonosatJNI.newSolver(args);
        } else {
            solverPtr = MonosatJNI.newSolver();
        }
        assert (solverPtr != 0);
        if (!enable_preprocessing) {
            disablePreprocessing();
        }
        true_lit = toLit(MonosatJNI.true_lit(solverPtr));
        initBV();
        initBuffers();
        Logic.setSolver(this);//default the current solver to the most recently created one.
    }

    private static String collectArgs(ArrayList<String> args) {
        String arg = "";
        for (String s : args) {
            arg += s + " ";
        }
        return arg;
    }

    @Override
    public synchronized void close() {
        //Does this method actually need to be syncronized?
        if (solverPtr != 0) {
            MonosatJNI.deleteSolver(solverPtr);
            solverPtr = 0;
        }
    }

    @Override
    protected void finalize() {
        //Consider replacing with Java 9's java.lang.ref.Cleaner in the future.
        //But for now, sticking with finalize (to support Java 8)
        close();
    }

    private void initBuffers() {
        initBuffer(0);
        initBuffer(1);
        initBuffer(2);
    }

    private void initBuffer(int bufferN) {
        if (bufferN == 0) {
            ByteBuffer b = ByteBuffer.allocateDirect(buffer_size0*4);//4 bytes per integer
            b.order(ByteOrder.LITTLE_ENDIAN);
            ints0 = b.asIntBuffer();
        } else if (bufferN == 1) {
            ByteBuffer b = ByteBuffer.allocateDirect(buffer_size0*4);
            b.order(ByteOrder.LITTLE_ENDIAN);
            ints1 = b.asIntBuffer();
        } else if (bufferN == 2) {
            ByteBuffer b = ByteBuffer.allocateDirect(buffer_size0*4);
            b.order(ByteOrder.LITTLE_ENDIAN);
            ints2 = b.asIntBuffer();
        }
    }

    protected IntBuffer getBuffer(int bufferN, int minsize) {

        if (bufferN == 0) {
            if (minsize >= buffer_size0) {
                buffer_size0 = minsize * 2;
                initBuffer(bufferN);
            }
            return ints0;
        } else if (bufferN == 1) {
            if (minsize >= buffer_size1) {
                buffer_size1 = minsize * 2;
                initBuffer(bufferN);
            }
            return ints1;
        } else if (bufferN == 2) {
            if (minsize >= buffer_size2) {
                buffer_size2 = minsize * 2;
                initBuffer(bufferN);
            }
            return ints2;
        }
        throw new RuntimeException("BufferN must be between 0 and 3");
    }

    /**
     * Instantiate a bitvector theory this
     */
    private void initBV() {
        assert (bvPtr == 0);
        bvPtr = MonosatJNI.initBVTheory(solverPtr);
        assert (bvPtr != 0);
    }

    public void setOutputFile(String file) {
        MonosatJNI.setOutputFile(solverPtr, file);
    }

    public Lit newLit() {
        return newLit(true);
    }

    public Lit newLit(boolean decisionVar) {
        int var = MonosatJNI.newVar(solverPtr);
        assert (var >= 0);
        Lit l = new Lit(var);
        while (lits.size() <= var) {
            lits.add(null);
        }
        assert (!l.sign());
        lits.set(l.toVar(), l);
        assert (toLit(var * 2) == l);
        return l;
    }

    void releaseLiteral(Lit l) {
        l.validate();
        int x = l.l;
        assert (x >= 0);
        MonosatJNI.releaseLiteral(solverPtr, x);
    }

    void setDecisionLiteral(Lit l, boolean decidable) {
        MonosatJNI.setDecisionVar(solverPtr, l.toVar(), decidable);
    }

    boolean isDecisionLiteral(Lit l) {
        return MonosatJNI.isDecisionVar(solverPtr, l.toVar());
    }

    public void setDecisionPriority(Lit l, int priority) {
        MonosatJNI.setDecisionPriority(solverPtr, l.toVar(), priority);
    }

    public int getDecisionPriority(Lit l) {
        return MonosatJNI.getDecisionPriority(solverPtr, l.toVar());
    }

    public void setDecisionPolarity(Lit l, boolean b) {
        MonosatJNI.setDecisionPolarity(solverPtr, l.toVar(), b);
    }

    public boolean getDecisionPolarity(Lit l) {
        return MonosatJNI.getDecisionPolarity(solverPtr, l.toVar());
    }

    public Lit True() {
        return true_lit;
    }

    public Lit False() {
        return true_lit.negate();
    }

    public void disallowSimplification(Lit l) {
        MonosatJNI.disallowLiteralSimplification(solverPtr, l.toVar());
    }

    public void disablePreprocessing() {
        MonosatJNI.disablePreprocessing(solverPtr);
    }

    public int nVars() {
        return MonosatJNI.nVars(solverPtr);
    }

    public int nClauses() {
        return MonosatJNI.nClauses(solverPtr);
    }

    public int nBitvectors() {
        return MonosatJNI.nBitvectors(solverPtr, bvPtr);
    }

    IntBuffer getBVBuffer(Collection<BitVector> clause, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, clause.size());
        int index = 0;
        for (BitVector bv : clause) {

            buffer.put(index, bv.id);
            index++;
        }
        return buffer;
    }

    IntBuffer getVarBuffer(Collection<Lit> clause, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, clause.size());
        int index = 0;
        for (Lit l : clause) {
            l.validate();
            buffer.put(index, l.toVar());
            index++;
        }
        return buffer;
    }

    IntBuffer getLitBuffer(Collection<Lit> clause) {
        return getLitBuffer(clause, 0);
    }

    IntBuffer getLitBuffer(Lit[] clause, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, clause.length);
        int index = 0;
        for (Lit l : clause) {
            l.validate();
            buffer.put(index, l.l);
            index++;
        }
        return buffer;
    }

    IntBuffer getLitBuffer(Collection<Lit> clause, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, clause.size());
        int index = 0;
        for (Lit l : clause) {
            l.validate();
            buffer.put(index, l.l);
            index++;
        }
        return buffer;
    }

    IntBuffer getIntBuffer(Collection<Integer> ints, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, ints.size());
        int index = 0;
        for (Integer i : ints) {
            buffer.put(index, i);
            index++;
        }
        return buffer;
    }

    /**
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     *
     * @param clause The clause to add to the solver
     * @return
     */
    public boolean addClause(ArrayList<Lit> clause) {

        boolean status = MonosatJNI.addClause(solverPtr, getLitBuffer(clause, 1), clause.size());
        return status;
    }

    //Solver API

    /**
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     *
     * @param l The literal to add as a unit clause to the solver (forcing l to be true)
     * @return
     */
    public boolean addClause(Lit l) {
        l.validate();
        boolean status = MonosatJNI.addUnitClause(solverPtr, l.l);
        return status;
    }

    public boolean addClause(Lit l1, Lit l2) {
        l1.validate();
        l2.validate();
        boolean status = MonosatJNI.addBinaryClause(solverPtr, l1.l, l2.l);
        return status;
    }

    public boolean addClause(Lit l1, Lit l2, Lit l3) {
        l1.validate();
        l2.validate();
        l3.validate();
        boolean status = MonosatJNI.addTertiaryClause(solverPtr, l1.l, l2.l, l3.l);
        return status;
    }

    //basic this functions
    public boolean solve() {
        boolean r = MonosatJNI.solve(solverPtr);
        return r;
    }

    public boolean solve(Lit... assumptions) {
        return MonosatJNI.solveAssumptions(solverPtr, getLitBuffer(assumptions, 0), assumptions.length);
    }

    public boolean solve(Collection<Lit> assumptions) {
        return MonosatJNI.solveAssumptions(solverPtr, getLitBuffer(assumptions), assumptions.size());
    }

    //Sets the (approximate) time limit in seconds before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable time limit.
    public void setTimeLimit(int seconds) {
        MonosatJNI.setTimeLimit(solverPtr, seconds);
    }

    //Sets the (approximate) memory limit in megabytes before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable memory limit.
    public void setMemoryLimit(int mb) {
        MonosatJNI.setMemoryLimit(solverPtr, mb);
    }

    //Sets the maximum number of (additional) conflicts allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable conflict limit.
    public void setConflictLimit(int num_conflicts) {
        MonosatJNI.setConflictLimit(solverPtr, num_conflicts);
    }

    //Sets the maximum number of (additional) propagations allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable propagation limit.
    public void setPropagationLimit(int num_propagations) {
        MonosatJNI.setPropagationLimit(solverPtr, num_propagations);
    }

    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    //Consider using Optional<Boolean> instead.
    public LBool solveLimited() {
        int result = MonosatJNI.solveLimited(solverPtr);
        assert (result >= 0);
        assert (result <= 2);
        return LBool.toLbool(result);
    }

    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    public LBool solveAssumptionsLimited(ArrayList<Lit> assumptions) {
        int result = MonosatJNI.solveAssumptionsLimited(solverPtr, getLitBuffer(assumptions), assumptions.size());
        assert (result >= 0);
        assert (result <= 2);
        return LBool.toLbool(result);
    }

    public boolean lastSolutionWasOptimal() {
        return MonosatJNI.lastSolutionWasOptimal(solverPtr);
    }

    /**
     * If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
     */
    public ArrayList<Lit> getConflictClause(){
        return getConflictClause(false);
    }

    /**
     * If the last solve call was UNSAT due to one or more assumption literals,
     * this method returns a subset of those assumption literals that are sufficient to keep make the instance UNSAT.
     * @param minimize If true, a (possibly expensive) search will be applied in the solver to find a locally minimal set of
     * literals that are mutually UNSAT. Else, an inexpensive, best effort set of literals will be returned.
     * @return
     */
    public ArrayList<Lit> getConflictClause(boolean minimize) {
        ArrayList<Lit> store = new ArrayList<Lit>();
        if(minimize){
            minimizeConflictClause();
        }

        //If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
        //Fills the given pointer with the first max_store_size literals of the conflict clause, and returns the number of literals in the conflict clause. Set store_clause to null and max_store_size to 0 to find the size of the conflict clause
        //Returns -1 if the solver has no conflict clause from the most recent solve() call (because that call was not UNSAT)
        int conflict_size = MonosatJNI.getConflictClause(solverPtr, null, 0);
        IntBuffer buf = getBuffer(0, conflict_size);
        int sz = MonosatJNI.getConflictClause(solverPtr, buf, conflict_size);
        assert (sz == conflict_size);
        for (int i = 0; i < conflict_size; i++) {
            int l = buf.get(i);
            Lit lit = toLit(l);
            store.add(lit);
        }
        return store;
    }


    //Backtrack the solver to level 0
    void backtrack() {
        MonosatJNI.backtrack(solverPtr);
    }

    protected Lit toLit(int literal) {
        assert (literal >= 0);
        int var = literal / 2;
        assert (var < nVars());//the variable must have already been declared in the sat solver before this call
        while (var >= lits.size()) {
            lits.add(null);
        }
        if (lits.get(var) == null) {
            lits.set(var, new Lit(var));
        }
        Lit l = lits.get(var);
        if ((literal & 1) == 1) {
            return l.negate();
        } else {
            return l;
        }
    }


    //Optimization API


    /**
     * Given a set of assumptions which are mutualy UNSAT, find a locally minimal subset that remains UNSAT.
     * (leaves the original set intact if the literals are not mutualy UNSAT)
     */
    public ArrayList<Lit> minimizeUnsatCore(Lit... literals){
        return minimizeUnsatCore(Arrays.asList(literals));
    }
    /**
     * Given a set of assumptions which are mutualy UNSAT, find a locally minimal subset that remains UNSAT.
     * (leaves the original set intact if the literals are not mutualy UNSAT)
     */
    public ArrayList<Lit> minimizeUnsatCore(Collection<Lit> assumptions) {
        ArrayList<Lit> store = new ArrayList<Lit>();
        IntBuffer buf = getLitBuffer(assumptions);
        int core_size = MonosatJNI.minimizeUnsatCore(solverPtr,buf,assumptions.size());
        assert(core_size>=0);
        assert(core_size<=assumptions.size());
      //  assumptions.clear();
        for (int i = 0; i < core_size; i++) {
            int l = buf.get(i);
            Lit lit = toLit(l);
            store.add(lit);
        }
        return store;
    }



    /**
     * After UNSAT solve calls with assumptions, the solver will find a 'conflict clause' consisting of a subset of the assumptions
     //which are sufficient to make the solver UNSAT (see getConflictClause).
     Normally, the conflict clause is produced as a side effect of proving the query unsat, with the solver only removing literals
     from the conflict clause on a best effort basis.
     This method will make repeated (and potentially expensive) calls to the SAT solver to attempt to remove further literals from
     the conflict clause.
     Afterward, the conflict clause can be obtained using getConflictClause().
     NOTE: this function may be expensive, is not required to get a conflict clause; getConflictClause() can be used after any UNSAT call with assumptions, even
     without calling minimizeConflictClause().
     Note also that if any of the setXXXXXLimits() are applied to the solver, this may not produce a locally minimal conflict clause.
     */
    private void minimizeConflictClause(){
        MonosatJNI.minimizeConflictClause(solverPtr);
    }

    public void clearOptimizationObjectives() {
        MonosatJNI.clearOptimizationObjectives(solverPtr);
    }

    public void maximizeBV(BitVector bv) {
        MonosatJNI.maximizeBV(solverPtr, bvPtr, bv.id);
    }

    public void minimizeBV(BitVector bv) {
        MonosatJNI.minimizeBV(solverPtr, bvPtr, bv.id);
    }

    public void maximizeLits(Collection<Lit> literals) {
        MonosatJNI.maximizeLits(solverPtr, getLitBuffer(literals), literals.size());
    }

    public void minimizeLits(Collection<Lit> literals) {
        MonosatJNI.minimizeLits(solverPtr, getLitBuffer(literals), literals.size());
    }

    public void maximizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        assert (literals.size() == weights.size());
        MonosatJNI.maximizeWeightedLits(solverPtr, getLitBuffer(literals), getIntBuffer(weights, 1), literals.size());
    }

    public void minimizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        assert (literals.size() == weights.size());
        MonosatJNI.minimizeWeightedLits(solverPtr, getLitBuffer(literals), getIntBuffer(weights, 1), literals.size());
    }

    public void assertAtMostOne(Collection<Lit> clause) {
        MonosatJNI.at_most_one(solverPtr, getLitBuffer(clause), clause.size());
    }

    public void assertPB(Collection<Lit> clause,  Comparison c, int compareTo) {
        assertPB(clause, null,c, compareTo);
    }

    public void assertPB(Collection<Lit> clause, Collection<Integer> weights,Comparison c, int compareTo) {
        IntBuffer wt_buffer = getBuffer(1, clause.size());
        int n_wts = 0;
        if (weights != null) {
            n_wts = Math.min(clause.size(), weights.size());
            int i = 0;
            for(Integer w:weights){
                wt_buffer.put(i, w);
                if(++i>=weights.size()){
                    break;
                }
            }
        }
        for (int i = n_wts; i < clause.size(); i++) {
            wt_buffer.put(i, 1); //default weight is 1
        }
        switch (c) {
            case LT:
                MonosatJNI.assertPB_lt(solverPtr, compareTo, clause.size(), getLitBuffer(clause), wt_buffer);
            case LEQ:
                MonosatJNI.assertPB_leq(solverPtr, compareTo, clause.size(), getLitBuffer(clause), wt_buffer);
            case EQ:
                MonosatJNI.assertPB_eq(solverPtr, compareTo, clause.size(), getLitBuffer(clause), wt_buffer);
            case GEQ:
                MonosatJNI.assertPB_geq(solverPtr, compareTo, clause.size(), getLitBuffer(clause), wt_buffer);
            case GT:
                MonosatJNI.assertPB_gt(solverPtr, compareTo, clause.size(), getLitBuffer(clause), wt_buffer);
        }
    }

    //Convert any pb constraints in the solver into cnf (will be called automatically before solve())
    public void flushPB() {
        MonosatJNI.flushPB(solverPtr);
    }

    /**
     * Caches constant bitvectors of value < MAX_CACHE_CONST
     *
     * @param width
     * @param constant
     * @return
     */
    public BitVector bv(int width, long constant) {
        if (constant >= 0 && constant <= MAX_CACHE_CONST) {
            while (cached_bvs.size() <= width) {
                cached_bvs.add(new ArrayList<>());
            }

            int small_const = (int) constant;
            while (cached_bvs.get(width).size() <= constant) {
                cached_bvs.get(width).add(null);
            }
            if (cached_bvs.get(width).get(small_const) == null) {
                cached_bvs.get(width).set(small_const, new BitVector(this, width, constant));
            }
            assert (cached_bvs.get(width).get(small_const) != null);
            return cached_bvs.get(width).get(small_const);
        }
        return new BitVector(this, width, constant);
    }

    public BitVector bv(int width) {
        return new BitVector(this, width);
    }

    /**
     * Clear the current satisfying model in the solver (if any), while preserving the current constraints.
     * There is normally no need to manually call restart(), as this is called automatically before calling solve().
     */
    public void restart(){
        MonosatJNI.backtrack(solverPtr);
    }

    /**
     * True if the solver has a satisfying model to its constraints
     * @return
     */
    public boolean hasModel(){
        return MonosatJNI.hasModel(solverPtr);
    }

    /**
     * Query the model in the solver, throwing an exception if the literal is unassigned in the model.
     * This can happen only if the literal is not a decision literal.
     */
    public boolean getValue(Lit l) throws RuntimeException {
        return getValue(l, LBool.Undef);
    }

    /**
     * Query the model in the solver.
     * If defaultVal is LBool.Undef, this will throw an exception if the literal is unassigned.
     * Else, if the literal is unassigned, defaultVal will be returned.
     */
    public boolean getValue(Lit l, LBool defaultVal) throws RuntimeException {
        LBool val = getPossibleValue(l);
        if (val == LBool.Undef) {
            if (defaultVal == LBool.Undef) {
                throw new RuntimeException("Literal " + l.toString() + " is unassigned in the current model");
            } else {
                val = defaultVal;
            }
        }
        return val == LBool.True;
    }

    /**
     * After a solve call, non-decision literals may or may not be assigned to a value.
     * Unassigned literals will have the value LBool.Undef;
     */
    public LBool getPossibleValue(Lit l) {
        return getPossibleValue(l,LBool.Undef);
    }

    /**
     * After a solve call, non-decision literals may or may not be assigned to a value.
     * Unassigned literals will have the value defaultValue
     */
    public LBool getPossibleValue(Lit l, LBool defaultValue) {
        l.validate();
        /*if (!MonosatJNI.hasModel(solverPtr)) {
            throw new RuntimeException("Solver has no model (this may indicate either that the solve() has not yet been called, or that the most recent call to solve() returned a value other than true, or that a constraint was added into the solver after the last call to solve()).");
        }*/
        LBool val = LBool.toLbool(MonosatJNI.getModel_Literal(solverPtr, l.l));
        if (val==LBool.Undef){
            return defaultValue;
        }else{
            return val;
        }
    }
    /**
     * Return the value of this bitvector from the solver.
     * Sometimes, a range of values may be determined by the solver to be satisfying.
     * If getMaximumValue is true, then largest value in that range will be returned,
     * otherwise, the smallest value is returned (this is relevant if optimization queries are being performed).
     *
     * @param bv
     * @param getMaximumValue
     * @return
     */
    public long getValue(BitVector bv, boolean getMaximumValue) {
        if (!MonosatJNI.hasModel(solverPtr)) {
            throw new RuntimeException("Solver has no model (this may indicate either that the solve() has not yet been called, or that the most recent call to solve() returned a value other than true, or that a constraint was added into the solver after the last call to solve()).");
        }
        return MonosatJNI.getModel_BV(solverPtr, bvPtr, bv.id, getMaximumValue);
    }

    /**
     * Return the value of this bitvector from the solver.
     * Sometimes, a range of values may be determined by the solver to be satisfying.
     * In this case, the smallest value is returned (this is relevant if optimization queries are being performed).
     *
     * @param bv
     * @return
     */
    public long getValue(BitVector bv) {
        return getValue(bv, false);
    }

    public Lit ite(Lit condition, Lit then, Lit els) {

        return this.toLit(MonosatJNI.Ite(this.solverPtr, condition.l, then.l, els.l));
    }

    //Higher level constructs for the solver
    //Literal level constructs

    public Lit and(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.And(this.solverPtr, args[0].l, args[1].l));
        }
        return and(Arrays.asList(args));
    }

    public Lit or(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Or(this.solverPtr, args[0].l, args[1].l));
        }
        return or(Arrays.asList(args));
    }

    public Lit not(Lit a) {
        return a.negate();
    }

    public Lit nand(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            return args[0].negate();
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Nand(this.solverPtr, args[0].l, args[1].l));
        }
        return nand(Arrays.asList(args));
    }

    public Lit nor(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Nor(this.solverPtr, args[0].l, args[1].l));
        }
        return nor(Arrays.asList(args));
    }

    public Lit xor(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            return args[0].negate();
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Xor(this.solverPtr, args[0].l, args[1].l));
        }
        return xor(Arrays.asList(args));
    }

    public Lit xnor(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Xnor(this.solverPtr, args[0].l, args[1].l));
        }
        return xnor(Arrays.asList(args));
    }

    //Assertion forms.
    public void assertTrue(Lit a) {
        MonosatJNI.Assert(this.solverPtr, a.l);
    }

    public void assertFalse(Lit a) {
        MonosatJNI.Assert(this.solverPtr, a.negate().l);
    }

    public void assertAnd(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            assertTrue(args[0]);
        } else if (args.length == 2) {
            MonosatJNI.AssertAnd(this.solverPtr, args[0].l, args[1].l);
        }
        assertAnd(Arrays.asList(args));
    }

    public void assertOr(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            assertTrue(args[0]);
        } else if (args.length == 2) {
            MonosatJNI.AssertOr(this.solverPtr, args[0].l, args[1].l);
        }
        assertOr(Arrays.asList(args));
    }

    public void assertNand(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            assertTrue(args[0].negate());
        } else if (args.length == 2) {
            MonosatJNI.AssertNand(this.solverPtr, args[0].l, args[1].l);
        }
        assertNand(Arrays.asList(args));
    }

    public void assertNor(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            assertTrue(args[0].negate());
        } else if (args.length == 2) {
            MonosatJNI.AssertNor(this.solverPtr, args[0].l, args[1].l);
        }
        assertNor(Arrays.asList(args));
    }

    public void assertXor(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            assertTrue(args[0]);//If this the correct behaviour here?
        } else if (args.length == 2) {
            MonosatJNI.AssertXor(this.solverPtr, args[0].l, args[1].l);
        }
        assertXor(Arrays.asList(args));
    }

    public void assertXnor(Lit... args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        } else if (args.length == 1) {
            assertTrue(args[0]);//If this the correct behaviour here?
        } else if (args.length == 2) {
            MonosatJNI.AssertXnor(this.solverPtr, args[0].l, args[1].l);
        }
        assertXnor(Arrays.asList(args));
    }

    public void assertEqual(Lit a, Lit b) {
        assertXnor(a, b);
    }
    public void assertImplies(Lit a, Lit b) {
        assertOr(a.negate(), b);
    }
    //Multi-Literal constructs
    public Lit and(Collection<Lit> elements) {

        return this.toLit(MonosatJNI.Ands(this.solverPtr, this.getLitBuffer(elements), elements.size()));
    }

    public Lit or(Collection<Lit> elements) {

        return this.toLit(MonosatJNI.Ors(this.solverPtr, this.getLitBuffer(elements), elements.size()));
    }

    public Lit nand(Collection<Lit> elements) {

        return this.toLit(MonosatJNI.Nands(this.solverPtr, this.getLitBuffer(elements), elements.size()));
    }

    public Lit nor(Collection<Lit> elements) {

        return this.toLit(MonosatJNI.Nors(this.solverPtr, this.getLitBuffer(elements), elements.size()));
    }

    public Lit xor(Collection<Lit> elements) {

        return this.toLit(MonosatJNI.Xors(this.solverPtr, this.getLitBuffer(elements), elements.size()));
    }

    public Lit xnor(Collection<Lit> elements) {

        return this.toLit(MonosatJNI.Xnors(this.solverPtr, this.getLitBuffer(elements), elements.size()));
    }

    //assertion forms
    public void assertAnd(Collection<Lit> elements) {

        MonosatJNI.AssertAnds(this.solverPtr, this.getLitBuffer(elements), elements.size());
    }

    public void assertOr(Collection<Lit> elements) {

        MonosatJNI.AssertOrs(this.solverPtr, this.getLitBuffer(elements), elements.size());
    }

    public void assertNand(Collection<Lit> elements) {

        MonosatJNI.AssertNands(this.solverPtr, this.getLitBuffer(elements), elements.size());
    }

    public void assertNor(Collection<Lit> elements) {

        MonosatJNI.AssertNors(this.solverPtr, this.getLitBuffer(elements), elements.size());
    }

    public void assertXor(Collection<Lit> elements) {

        MonosatJNI.AssertXors(this.solverPtr, this.getLitBuffer(elements), elements.size());
    }

    public void assertXnor(Collection<Lit> elements) {

        MonosatJNI.AssertXnors(this.solverPtr, this.getLitBuffer(elements), elements.size());
    }

    //Bitvector constructs
    public BitVector ite(Lit condition, BitVector then, BitVector els) {

        assert (then.width() == els.width());
        BitVector result = new BitVector(this, then.width());
        MonosatJNI.bv_ite(this.solverPtr, this.bvPtr, condition.l, then.id, els.id, result.id);
        return result;
    }

    public BitVector and(BitVector a, BitVector b) {

        return a.and(b);
    }

    public BitVector or(BitVector a, BitVector b) {
        return a.or(b);
    }

    public BitVector not(BitVector a) {
        return a.not();
    }

    public BitVector nand(BitVector a, BitVector b) {
        return a.nand(b);
    }

    public BitVector nor(BitVector a, BitVector b) {
        return a.nor(b);
    }

    public BitVector xor(BitVector a, BitVector b) {
        return a.xor(b);
    }

    public BitVector xnor(BitVector a, BitVector b) {
        return a.xnor(b);
    }

    public BitVector add(BitVector a, BitVector b) {
        return a.add(b);
    }

    public BitVector subtract(BitVector a, BitVector b) {
        return a.subtract(b);
    }

    public void assertEqual(BitVector a, BitVector b) {
        assertTrue(a.geq(b));
        assertTrue(a.leq(b));
    }

    public void assertEqual(BitVector a, long constant) {
        BitVector b = this.bv(a.width(), constant);
        assertTrue(a.geq(b));
        assertTrue(a.leq(b));
    }

    public void assertEqual(long constant, BitVector a) {
        BitVector b = this.bv(a.width(), constant);
        assertTrue(a.geq(b));
        assertTrue(a.leq(b));
    }

    public BitVector min(Collection<BitVector> args) {

        assert (args.size() >= 0);
        int w = args.iterator().next().width();
        BitVector result = new BitVector(this, w);
        MonosatJNI.bv_min(this.solverPtr, this.bvPtr, this.getBVBuffer(args, 0), args.size(), result.id);
        return result;
    }

    public BitVector min(BitVector a, BitVector b) {
        ArrayList pair = new ArrayList();
        pair.add(a);
        pair.add(b);
        return min(pair);
    }

    public BitVector max(Collection<BitVector> args) {

        assert (args.size() >= 0);
        int w = args.iterator().next().width();
        BitVector result = new BitVector(this, w);
        MonosatJNI.bv_min(this.solverPtr, this.bvPtr, this.getBVBuffer(args, 0), args.size(), result.id);
        return result;
    }

    public BitVector max(BitVector a, BitVector b) {
        ArrayList pair = new ArrayList();
        pair.add(a);
        pair.add(b);
        return max(pair);
    }




}
