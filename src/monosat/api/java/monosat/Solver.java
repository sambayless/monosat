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


import java.io.Closeable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.*;

import java.util.logging.Logger;

/**
 * Represents a MonoSAT solver instance.
 * Multiple solvers may be instantiated, with separate
 * sets of literals and clauses.
 */
public final class Solver implements Closeable {
    /**
     * Holds weak references to all currently existing solvers,
     * so that global logic operations on True/False can be applied.
     */
    protected static WeakHashMap<Solver,Boolean> solvers = new WeakHashMap<Solver, Boolean>();
    protected static final Logger log = Logger.getLogger("monosat");

    /**
     * Largest constant BitVector value to cache.
     * This may change in the future.
     */
    private static long MAX_CACHE_CONST = 255;

    /**
     * Handle to the underlying monosat solver instance.
     * This is really a pointer, masquerading as a long.
     */
    protected long solverPtr = 0;

    /**
     * Handle to the underlying monosat BitVector theory instance.
     * This is really a pointer, masquerading as a long.
     */
    protected long bvPtr = 0;

    //Used internally to manage byte buffers for calls to the C api
    private int buffer_size0 = 1024;
    private int buffer_size1 = 1024;
    private int buffer_size2 = 1024;

    private IntBuffer ints0;
    private IntBuffer ints1;
    private IntBuffer ints2;

    //Caches instantiated, small BitVectors
    private ArrayList<ArrayList<BitVector>> cached_bvs = new ArrayList<ArrayList<BitVector>>();

    //Holds instances of all literals, so that we don't need to create multiple literal objects for the same literal
    private ArrayList<Lit> allLits = new ArrayList<>();

    //contains only the positive versions of instianted literals, in the order they were created.
    private LinkedHashSet<Lit> positiveLiterals = new LinkedHashSet<>();

    private ArrayList<BitVector> allBVs = new ArrayList<>();

    /**
     * Represents a value that is either true, false, or undefined.
     */
    protected enum LBool {
        //Don't change the order of these, as they must match the order of the l_bool enum defined in Monosat.
        True, False, Undef;

        public static LBool toLbool(int value) {
            return values()[value];
        }
        public static LBool fromBool(boolean value) {
            return values()[value?0:1];
        }

        /**
         * Optional value objects are created statically allocated for each LBool type,
         * so that we don't need to instantiate new Optional instances on the heap for toOpt() calls.
         * This matters, as a user may look up the values of millions of literals in common use cases.
         */
        @SuppressWarnings("OptionalUsedAsFieldOrParameterType")
        private Optional<Boolean> opt;

        static {
            True.opt = Optional.of(true);
            False.opt = Optional.of(false);
            Undef.opt = Optional.empty();
        }

        public Optional<Boolean> toOpt(){
            return opt;
        }

    }


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

    public Solver(boolean enablePreprocessing) {
        this("", enablePreprocessing);
    }

    public Solver(ArrayList<String> args) {
        this(collectArgs(args), false);
    }

    public Solver(ArrayList<String> args, String outputFile) {
        this(collectArgs(args), false,outputFile);
    }

    public Solver(ArrayList<String> args, boolean enablePreprocessing) {
        this(collectArgs(args), enablePreprocessing,"");
    }
    public Solver(ArrayList<String> args, boolean enablePreprocessing, String outputFile) {
        this(collectArgs(args), enablePreprocessing,outputFile);
    }
    public Solver(String args, boolean enablePreprocessing){
        this(args,enablePreprocessing,"");
    }
    public Solver(String args, String outputFule){
        this(args,false,outputFule);
    }
    public Solver(String args, boolean enablePreprocessing,String outputFile) {
        if (args != null && args.length() > 0) {
            solverPtr = MonosatJNI.newSolver("monosat " + args);
        } else {
            solverPtr = MonosatJNI.newSolver();
        }
        if (solverPtr==0){
            throw new RuntimeException("Failed to created solver");
        }
        //Keep a global list of all solvers that yet created
        solvers.put(this,true);
        if(outputFile!=null && outputFile.length()>0){
            MonosatJNI.setOutputFile(solverPtr, outputFile);
        }
        if (!enablePreprocessing) {
            disablePreprocessing();
        }
        registerLit(Lit.True,Lit.False);
        this.addClause(Lit.True);
        //Ensure that the True lit returned by circuit operations (eg, and())
        //is the same as this True Lit.
        assert(Lit.True== toLit(MonosatJNI.true_lit(solverPtr)));

        initBV();
        initBuffers();
        assert(positiveLiterals.contains(Lit.True));
    }

    /**
     * @return The version string of the MonoSAT native library.
     */
    public static String getVersion(){
        return MonosatJNI.getVersion();
    }

    private static String collectArgs(ArrayList<String> args) {
        StringBuilder arg = new StringBuilder();
        for (String s : args) {
            arg.append(s).append(" ");
        }
        return arg.toString();
    }

    /**
     * Release all native resources associated with this solver.
     * This does not normally need to be called manually (though it is safe to do so),
     * as it will be called automatically during garbage collection.
     */
    @Override
    public synchronized void close() {
        //Does this method actually need to be syncronized?
        if (solverPtr != 0) {
            solvers.remove(this);
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

    protected void validate(Lit l){
        if(l==Lit.True || l==Lit.False)
            return;
        if(l==null){
            throw new IllegalArgumentException("Literal is null");
        }else if (l.l<0 ) {
            throw new IllegalArgumentException("Literal " + l.toString() + " is not a valid literal.");
        }else if (l.solver!=this){
            throw new IllegalArgumentException("Cannot pass literal belonging to solver " + (l.solver ==null?"null":l.solver.toString()) + " to solver " + toString());
        }else if(l.toVar()>=nVars()){
            throw new IllegalArgumentException("Literal is undefined in solver (too large)");
        }
    }

    protected void validate(Lit... args){
        for(Lit l:args){
            validate(l);
        }
    }
    protected void validate(Collection<Lit> args){
        for(Lit l:args){
            validate(l);
        }
    }

    protected void validate(BitVector... args){
        for(BitVector bv:args){
            if(bv==null){
                throw new IllegalArgumentException("Literal is null");
            }else if (bv.getSolver() !=this){
                throw new IllegalArgumentException("Cannot pass literal belonging to solver " +  (bv.getSolver() == null ? "null" : bv.getSolver().toString()) + " to solver " + toString());
            }else if (bv.id<0 ) {
                throw new IllegalArgumentException("Bitvector is undefined " + bv.toString());
            }
        }
    }
    protected void validateBV(Collection<BitVector> args){
        for(BitVector bv:args){
            if(bv==null){
                throw new IllegalArgumentException("Literal is null");
            }else if (bv.getSolver() !=this){
                throw new IllegalArgumentException("Cannot pass literal belonging to solver " +  (bv.getSolver() == null ? "null" : bv.getSolver().toString()) + " to solver " + toString());
            }else if (bv.id<0 ) {
                throw new IllegalArgumentException("Bitvector is undefined " + bv.toString());
            }
        }
    }

    public void releaseLiteral(Lit l) {
        validate(l);
        int x = l.l;
        assert (x >= 0);
        MonosatJNI.releaseLiteral(solverPtr, x);
    }

    public void setDecisionLiteral(Lit l, boolean decidable) {
        validate(l);
        MonosatJNI.setDecisionVar(solverPtr, l.toVar(), decidable);
    }

    public boolean isDecisionLiteral(Lit l) {
        validate(l);
        return MonosatJNI.isDecisionVar(solverPtr, l.toVar());
    }

    public void setDecisionPriority(Lit l, int priority) {
        validate(l);
        MonosatJNI.setDecisionPriority(solverPtr, l.toVar(), priority);
    }

    public int getDecisionPriority(Lit l) {
        validate(l);
        return MonosatJNI.getDecisionPriority(solverPtr, l.toVar());
    }

    public void setDecisionPolarity(Lit l, boolean b) {
        validate(l);
        MonosatJNI.setDecisionPolarity(solverPtr, l.toVar(), b);
    }

    public boolean getDecisionPolarity(Lit l) {
        validate(l);
        return MonosatJNI.getDecisionPolarity(solverPtr, l.toVar());
    }


    public void disallowSimplification(Lit l) {
        validate(l);
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


    /**
     * Internal method for converting java arrays of bitvectors into
     * byte buffers (as required by the Monosat JNI)
     * @param clause
     * @param bufferN
     * @return
     */
    protected IntBuffer getBVBuffer(Collection<BitVector> clause, int bufferN) {
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

    /**
     * Internal method for converting java collections of allLits into
     * byte buffers (as required by the Monosat JNI)
     * @param clause
     * @param bufferN
     * @return
     */
    protected IntBuffer getVarBuffer(Collection<Lit> clause, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, clause.size());
        int index = 0;
        for (Lit l : clause) {
            buffer.put(index, l.toVar());
            index++;
        }
        return buffer;
    }

    /**
     * Internal method for converting java collections of Lits into
     * byte buffers (as required by the Monosat JNI)
     * @param clause
     * @return
     */
    protected IntBuffer getLitBuffer(Collection<Lit> clause) {
        return getLitBuffer(clause, 0);
    }

    /**
     * Internal method for converting java arrays of Lits into
     * byte buffers (as required by the Monosat JNI)
     * @param clause
     * @param bufferN
     * @return
     */
    protected IntBuffer getLitBuffer(Lit[] clause, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, clause.length);
        int index = 0;
        for (Lit l : clause) {
            buffer.put(index, l.l);
            index++;
        }
        return buffer;
    }

    /**
     * Internal method for converting java collections of Lits into
     * byte buffers (as required by the Monosat JNI)
     * @param clause
     * @param bufferN
     * @return
     */
    protected IntBuffer getLitBuffer(Collection<Lit> clause, int bufferN) {
        assert (bufferN < 3);
        assert (bufferN >= 0);
        IntBuffer buffer = getBuffer(bufferN, clause.size());
        int index = 0;
        for (Lit l : clause) {
            buffer.put(index, l.l);
            index++;
        }
        return buffer;
    }

    /**
     * Internal method for converting java collections of integers into
     * byte buffers (as required by the Monosat JNI)
     * @param ints
     * @param bufferN
     * @return
     */
    protected IntBuffer getIntBuffer(Collection<Integer> ints, int bufferN) {
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

    //Solver API

    /**
     * Add a clause to the solver.
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     */
    public boolean addClause(Lit a) {
        validate(a);
        return MonosatJNI.addUnitClause(solverPtr, a.l);
    }

    /**
     * Add a clause to the solver.
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     */
    public boolean addClause(Lit a, Lit b) {
        validate(a);
        validate(b);
        return MonosatJNI.addBinaryClause(solverPtr, a.l,b.l);
    }

    /**
     * Add a clause to the solver.
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     */
    public boolean addClause(Lit a, Lit b,Lit c) {
        validate(a);
        validate(b);
        validate(c);
        return MonosatJNI.addTertiaryClause(solverPtr, a.l,b.l,c.l);
    }

    /**
     * Add a clause to the solver.
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     *
     * @param args The literals to add as a clause to the solver (forcing at least one of them to be true)
     * @return
     */
    public boolean addClause(Lit... args) {
        validate(args);
        return MonosatJNI.addClause(solverPtr, getLitBuffer(args, 1), args.length);
    }

    /**
     * Add a clause to the solver.
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     *
     * @param clause The clause to add to the solver
     * @return
     */
    public boolean addClause(Collection<Lit> clause) {
        validate(clause);
        return MonosatJNI.addClause(solverPtr, getLitBuffer(clause, 1), clause.size());
    }

    //basic this functions

    /**
     * Either find a satisfying solution to the constraints in te formula, or
     * prove the formula to be unsatisfiable.
     * @return True if a satisfying soltuion was found, false if it does not.
     */
    public boolean solve() {
        return MonosatJNI.solve(solverPtr);
    }
    /**
     * Either find a satisfying solution to the constraints in te formula, or
     * prove the formula to be unsatisfiable, while temporarily enforcing the
     * literals in assumptions to be true.     * *
     * @return True if a satisfying soltuion was found, false if it does not.
     */
    public boolean solve(Lit... assumptions) {
        validate(assumptions);
        return MonosatJNI.solveAssumptions(solverPtr, getLitBuffer(assumptions, 0), assumptions.length);
    }
    /**
     * Either find a satisfying solution to the constraints in te formula, or
     * prove the formula to be unsatisfiable, while temporarily enforcing the
     * literals in assumptions to be true.     * *
     * @return True if a satisfying soltuion was found, false if it does not.
     */
    public boolean solve(Collection<Lit> assumptions) {
        validate(assumptions);
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
    /**
     * Attempt to find a satisfying solution, but return Optional.empty() if any resource limits
     * are violated. To set resource limits, see:
     * setTimeLimit()
     * setMemoryLimit()
     * setConflictLimit()
     * setPropagationLimit()
    */
    public Optional<Boolean> solveLimited() {
        int result = MonosatJNI.solveLimited(solverPtr);
        assert (result >= 0);
        assert (result <= 2);
        return LBool.toLbool(result).toOpt();
    }

    /**
     * Attempt to find a satisfying solution, but return Optional.empty() if any resource limits
     * are violated. To set resource limits, see:
     * setTimeLimit()
     * setMemoryLimit()
     * setConflictLimit()
     * setPropagationLimit()
     */
    public Optional<Boolean> solveLimited(Collection<Lit> assumptions) {
        validate(assumptions);
        int result = MonosatJNI.solveAssumptionsLimited(solverPtr, getLitBuffer(assumptions), assumptions.size());
        assert (result >= 0);
        assert (result <= 2);
        return LBool.toLbool(result).toOpt();
    }
    /**
     * Attempt to find a satisfying solution, but return Optional.empty() if any resource limits
     * are violated. To set resource limits, see:
     * setTimeLimit()
     * setMemoryLimit()
     * setConflictLimit()
     * setPropagationLimit()
     */
    public Optional<Boolean> solveLimited(Lit... assumptions) {
        validate(assumptions);
        int result = MonosatJNI.solveAssumptionsLimited(solverPtr,  getLitBuffer(assumptions, 0), assumptions.length);
        assert (result >= 0);
        assert (result <= 2);
        return LBool.toLbool(result).toOpt();
    }
    /**
     * Returns true if the last solve() call found an optimal solution.
     * This is normally always the case, unless resource limits are enforced.
     * @return
     */
    public boolean lastSolutionWasOptimal() {
        return MonosatJNI.lastSolutionWasOptimal(solverPtr);
    }

    /**
     * If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
     */
    public List<Lit> getConflictClause(){
        return getConflictClause(false);
    }

    /**
     * If the last solve call was UNSAT due to one or more assumption literals,
     * this method returns a subset of those assumption literals that are sufficient to keep make the instance UNSAT.
     * @param minimize If true, a (possibly expensive) search will be applied in the solver to find a locally minimal set of
     * literals that are mutually UNSAT. Else, an inexpensive, best effort set of literals will be returned.
     * @return
     */
    public List<Lit> getConflictClause(boolean minimize) {
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
        validate(store);
        return store;
    }




    /**
     * Internal method used to cache newly instantiated literals,
     * so that literals returned from C++ can reuse the same objects.
     * @param l
     */
    protected void registerLit(Lit l){
        registerLit(l,null);
    }


    /**
     * Internal method used to maintain an iterable list of bitvectors in the solver.
     *
     */
    protected void registerBitVector(BitVector bv){
        allBVs.add(bv);
    }


    /**
     * Internal method used to cache newly instantiated literals,
     * so that literals returned from C++ can reuse the same objects.
     * @param l
     */
    protected void registerLit(Lit l, Lit notL){
        assert(l.l>=0);
        int literal = l.toInt();
        assert (literal >= 0);
        int var = literal / 2;
        assert (var < nVars());//the variable must have already been declared in the sat solver before this call
        while (var*2+1 >= allLits.size()) {
            allLits.add(null);
        }

        assert (allLits.get(var*2) == null);
        assert (allLits.get(var*2+1) == null);

        if(notL==null){
            if(l.sign()){
                notL = new Lit(this,var*2);
            }else{
                notL = new Lit(this,var*2+1);
            }
        }

        if(l.sign()) {
            assert(notL.l==l.l-1);
            allLits.set(var*2, notL);
            allLits.set(var*2+1,l);
        }else {
            assert(notL.l==l.l+1);
            allLits.set(var*2,l);
            allLits.set(var*2+1,notL);
        }

        assert(!allLits.get(var*2).sign());
        assert(allLits.get(var*2+1).sign());

        assert(allLits.get(l.l)!=null);
        assert(allLits.get(l.l)==l);

        if(l.sign()){
            assert (!positiveLiterals.contains(notL));
            positiveLiterals.add(notL);
        }else {
            assert (!positiveLiterals.contains(l));
            positiveLiterals.add(l);
        }
    }

    /**
     * Internal method to convert C++ literal integers into java Lits
     */
    protected Lit toLit(int literal) {
        assert (literal >= 0);
        int var = literal / 2;
        assert (var < nVars());//the variable must have already been declared in the sat solver before this call
        while (var*2+1 >= allLits.size()) {
            allLits.add(null);
        }
        if (allLits.get(var*2) == null) {
            assert(allLits.get(var*2+1) == null);
            Lit l =  new Lit(this,var*2);
            Lit notL = new Lit(this,var*2+1);
            allLits.set(var*2,l);
            allLits.set(var*2+1, notL);
            assert(!allLits.get(var*2).sign());
            assert(allLits.get(var*2+1).sign());

            assert (!positiveLiterals.contains(l));
            positiveLiterals.add(l);
        }
        assert(allLits.get(literal)!=null);
        assert(allLits.get(literal).l==literal);
        return allLits.get(literal);
    }

    /**
     * @return immutable view of the positive literals in the solver.
     */
    public Collection<Lit> getLits(){
        return Collections.unmodifiableSet(positiveLiterals);
    }

    /**
     * @return immutable view of the BitVectors in the solver.
     */
    public Collection<BitVector> getBitVectors(){
        return Collections.unmodifiableList(allBVs);
    }


    //Optimization API

    /**
     * Given a set of assumptions which are mutualy UNSAT, find a locally minimal subset that remains UNSAT.
     * (leaves the original set intact if the literals are not mutualy UNSAT)
     */
    public List<Lit> minimizeUnsatCore(Lit... literals){
        validate(literals);
        return minimizeUnsatCore(Arrays.asList(literals));
    }
    /**
     * Given a set of assumptions which are mutualy UNSAT, find a locally minimal subset that remains UNSAT.
     * (leaves the original set intact if the literals are not mutualy UNSAT)
     */
    public List<Lit> minimizeUnsatCore(Collection<Lit> literals) {
        validate(literals);
        ArrayList<Lit> store = new ArrayList<Lit>();
        IntBuffer buf = getLitBuffer(literals);
        int core_size = MonosatJNI.minimizeUnsatCore(solverPtr,buf,literals.size());
        assert(core_size>=0);
        assert(core_size<=literals.size());
        for (int i = 0; i < core_size; i++) {
            int l = buf.get(i);
            Lit lit = toLit(l);
            store.add(lit);
        }
        validate(store);
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

    /**
     * Clear any optimizaiton objectives in the solver.
     */
    public void clearOptimizationObjectives() {
        MonosatJNI.clearOptimizationObjectives(solverPtr);
    }

    /**
     * Add an optimization objective to the solver:
     * maximize the specified BitVector.
     */
    public void maximizeBV(BitVector bv) {
        validate(bv);
        MonosatJNI.maximizeBV(solverPtr, bvPtr, bv.id);
    }

    /**
     * Add an optimization objective to the solver:
     * minimize the specified BitVector.
     */
    public void minimizeBV(BitVector bv) {
        validate(bv);
        MonosatJNI.minimizeBV(solverPtr, bvPtr, bv.id);
    }

    /**
     * Add an optimization objective to the solver:
     * maximize the number of true literals from among the given literals.
     */
    public void maximizeLits(Collection<Lit> literals) {
        validate(literals);
        MonosatJNI.maximizeLits(solverPtr, getLitBuffer(literals), literals.size());
    }

    /**
     * Add an optimization objective to the solver:
     * minimize the number of true literals from among the given literals.
     */
    public void minimizeLits(Collection<Lit> literals) {
        validate(literals);
        MonosatJNI.minimizeLits(solverPtr, getLitBuffer(literals), literals.size());
    }

    /**
     * Add an optimization objective to the solver:
     * maximize the weighted number of true literals from among the given literals.
     */
    public void maximizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        validate(literals);
        if (literals.size()!=weights.size()){
            throw  new IllegalArgumentException("literals and weights must be the same length");
        }
        MonosatJNI.maximizeWeightedLits(solverPtr, getLitBuffer(literals), getIntBuffer(weights, 1), literals.size());
    }

    /**
     * Add an optimization objective to the solver:
     * minimize the weighted number of true literals from among the given literals.
     */
    public void minimizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        validate(literals);
        if (literals.size()!=weights.size()){
            throw  new IllegalArgumentException("literals and weights must be the same length");
        }
        MonosatJNI.minimizeWeightedLits(solverPtr, getLitBuffer(literals), getIntBuffer(weights, 1), literals.size());
    }

    /**
     * Enforce that at most one of the specified literals may be true.
     * If 1 or fewer arguments are given, has no effect.
     * If exactly 2 arguments are given, this is the same as:
     * assertOr(args[0],args[1])
     */
    public void assertAtMostOne(Lit... args) {
        validate(args);
        assertAtMostOne(Arrays.asList(args));
    }

    /**
     * Enforce that at most one of the specified literals may be true.
     * If 1 or fewer arguments are given, has no effect.
     * If exactly 2 arguments are given, this is the same as:
     * assertOr(args[0],args[1])
     */
    public void assertAtMostOne(Collection<Lit> args) {
        //simple at-most-one constraint: asserts that at most one of the set of variables (NOT LITERALS) may be true.
        //for small numbers of variables, consider using a direct CNF encoding instead
        validate(args);
        if(args.size()<=6){
            //make this constant configurable in the future
            //for small enough sets of literals, directly instantiate the binary clauses constraining them
            //rather than introduce an amo theory
            MonosatJNI.AssertAMO(solverPtr,getLitBuffer(args),args.size());
        }else{
            //workaround for internal limitations in monosat which require the
            //amo constraints to operate only on variables (not literals), which must also not have been used elsewhere
            ArrayList<Integer> vars = new ArrayList<Integer>();
            for(Lit l:args){
                Lit l2 = new Lit(this,false);
                assertEqual(l,l2);
                vars.add(l2.toVar());
            }
            MonosatJNI.at_most_one(solverPtr, getIntBuffer(vars,0), args.size());
        }
    }

    /**
     * Enforce a pseudo-Boolean constraint.
     * The number of true literals from among args must satisfy comparison c
     * relative to compare to.
     *
     * For example, if c is Comparison.LEQ, and compareTo is 3,
     * then at most 3 literals from args may be true.
     */
    public void assertPB(Collection<Lit> args,  Comparison c, int compareTo) {
        validate(args);
        ArrayList<Lit> tmp = new ArrayList<>();
        tmp.addAll(args);
        assertPB(tmp, null,c, compareTo);
    }


    /**
     * Enforce a weighted pseudo-Boolean constraint.
     * The weighted number of true literals from among args must satisfy comparison c
     * relative to compare to.
     *
     */
    public void assertPB(List<Lit> args, List<Integer> weights, Comparison c, int compareTo) {
        validate(args);

        IntBuffer wt_buffer = getBuffer(1, args.size());
        int n_wts = 0;
        if (weights != null) {
            if(weights.size()>args.size()){
                throw  new IllegalArgumentException("Must not have more weights then literals.");
            }
            n_wts = Math.min(args.size(), weights.size());
            int i = 0;
            for(Integer w:weights){
                wt_buffer.put(i, w);
                if(++i>=weights.size()){
                    break;
                }
            }
        }
        for (int i = n_wts; i < args.size(); i++) {
            wt_buffer.put(i, 1); //default weight is 1
        }
        switch (c) {
            case LT:
                MonosatJNI.assertPB_lt(solverPtr, compareTo, args.size(), getLitBuffer(args), wt_buffer);
            case LEQ:
                MonosatJNI.assertPB_leq(solverPtr, compareTo, args.size(), getLitBuffer(args), wt_buffer);
            case EQ:
                MonosatJNI.assertPB_eq(solverPtr, compareTo, args.size(), getLitBuffer(args), wt_buffer);
            case GEQ:
                MonosatJNI.assertPB_geq(solverPtr, compareTo, args.size(), getLitBuffer(args), wt_buffer);
            case GT:
                MonosatJNI.assertPB_gt(solverPtr, compareTo, args.size(), getLitBuffer(args), wt_buffer);
        }
    }
    /**
     * Immediately convert any pseudo Boolean constraints in the solver into clauses.
     * This function does not need to be manually called, as it will be called automatically before 'solve' calls.
    */
    protected void flushPB() {
        MonosatJNI.flushPB(solverPtr);
    }

    /**
     * Caches constant bitvectors of value <= MAX_CACHE_CONST
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

    /**
     * Reset any decisions in the solver, while preserving the current constraints.
     * There is normally no need to manually call restart(), as this is called automatically before calling solve().
     */
    public void restart(){
        MonosatJNI.backtrack(solverPtr);
    }

    /**
     * True if the solver has a satisfying model to its constraints.
     * @return
     */
    public boolean hasModel(){
        return MonosatJNI.hasModel(solverPtr);
    }




    public Lit ite(Lit condition, Lit then, Lit els) {
        validate(condition,then,els);
        return this.toLit(MonosatJNI.Ite(this.solverPtr, condition.l, then.l, els.l));
    }

    //Higher level constructs for the solver
    //Literal level constructs

    public Lit and(Lit... args) {
        validate(args);
        if (args.length == 0) {
            return Lit.True;
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.And(this.solverPtr, args[0].l, args[1].l));
        }
        return and(Arrays.asList(args));
    }

    public Lit or(Lit... args) {
        validate(args);
        if (args.length == 0) {
            return Lit.False;
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {
            return this.toLit(MonosatJNI.Or(this.solverPtr, args[0].l, args[1].l));
        }
        return or(Arrays.asList(args));
    }

    public Lit not(Lit a) {
        validate(a);
        int l = a.toInt();
        l = l^1;//bit twiddle odd to even
        assert(l>=0);
        assert(l< allLits.size());
        return allLits.get(l);
    }

    public Lit nand(Lit... args) {
        validate(args);
        if (args.length == 0) {
            return this.toLit(MonosatJNI.Nands(this.solverPtr, getBuffer(0,0),0));
        } else if (args.length == 1) {
            return args[0].not();
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Nand(this.solverPtr, args[0].l, args[1].l));
        }
        return nand(Arrays.asList(args));
    }

    public Lit nor(Lit... args) {
        validate(args);
        if (args.length == 0) {
            return this.toLit(MonosatJNI.Nors(this.solverPtr, getBuffer(0,0),0));
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {
            return this.toLit(MonosatJNI.Nor(this.solverPtr, args[0].l, args[1].l));
        }
        return nor(Arrays.asList(args));
    }

    public Lit xor(Lit... args) {
        validate(args);
        if (args.length == 0) {
            return this.toLit(MonosatJNI.Xors(this.solverPtr, getBuffer(0,0),0));
        } else if (args.length == 1) {
            return args[0].not();
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Xor(this.solverPtr, args[0].l, args[1].l));
        }
        return xor(Arrays.asList(args));
    }

    public Lit xnor(Lit... args) {
        validate(args);
        if (args.length == 0) {
            return this.toLit(MonosatJNI.Xnors(this.solverPtr, getBuffer(0,0),0));
        } else if (args.length == 1) {
            return args[0];
        } else if (args.length == 2) {

            return this.toLit(MonosatJNI.Xnor(this.solverPtr, args[0].l, args[1].l));
        }
        return xnor(Arrays.asList(args));
    }

    //Assertion forms.
    public void assertTrue(Lit a) {
        validate(a);
        MonosatJNI.Assert(this.solverPtr, a.l);
    }

    public void assertFalse(Lit a) {
        validate(a);
        MonosatJNI.Assert(this.solverPtr, a.not().l);
    }

    public void assertAnd(Lit... args) {
        validate(args);
        if (args.length == 0) {
            //do nothing
        } else if (args.length == 1) {
            assertTrue(args[0]);
        } else if (args.length == 2) {
            MonosatJNI.AssertAnd(this.solverPtr, args[0].l, args[1].l);
        }
        assertAnd(Arrays.asList(args));
    }

    public void assertOr(Lit... args) {
        validate(args);
        if (args.length == 0) {
            //this is a contradiction
            MonosatJNI.AssertOrs(this.solverPtr, getBuffer(0,0), 0);
        } else if (args.length == 1) {
            assertTrue(args[0]);
        } else if (args.length == 2) {
            MonosatJNI.AssertOr(this.solverPtr, args[0].l, args[1].l);
        }
        assertOr(Arrays.asList(args));
    }

    public void assertNand(Lit... args) {
        validate(args);
        if (args.length == 0) {
            MonosatJNI.AssertNands(this.solverPtr, getBuffer(0,0), 0);
        } else if (args.length == 1) {
            assertTrue(args[0].not());
        } else if (args.length == 2) {
            MonosatJNI.AssertNand(this.solverPtr, args[0].l, args[1].l);
        }
        assertNand(Arrays.asList(args));
    }

    public void assertNor(Lit... args) {
        validate(args);
        if (args.length == 0) {
            MonosatJNI.AssertNors(this.solverPtr, getBuffer(0,0), 0);
        } else if (args.length == 1) {
            assertTrue(args[0].not());
        } else if (args.length == 2) {
            MonosatJNI.AssertNor(this.solverPtr, args[0].l, args[1].l);
        }
        assertNor(Arrays.asList(args));
    }

    public void assertXor(Lit... args) {
        validate(args);
        if (args.length == 0) {
            MonosatJNI.AssertXors(this.solverPtr, getBuffer(0,0), 0);
        } else if (args.length == 1) {
            assertTrue(args[0]);//If this the correct behaviour here?
        } else if (args.length == 2) {
            MonosatJNI.AssertXor(this.solverPtr, args[0].l, args[1].l);
        }
        assertXor(Arrays.asList(args));
    }

    public void assertXnor(Lit... args) {
        validate(args);
        if (args.length == 0) {
            MonosatJNI.AssertXnors(this.solverPtr, getBuffer(0,0), 0);
        } else if (args.length == 1) {
            assertTrue(args[0]);//If this the correct behaviour here?
        } else if (args.length == 2) {
            MonosatJNI.AssertXnor(this.solverPtr, args[0].l, args[1].l);
        }
        assertXnor(Arrays.asList(args));
    }

    public void assertEqual(Lit a, Lit b) {
        validate(a,b);
        assertXnor(a, b);
    }
    public void assertImplies(Lit a, Lit b) {
        validate(a,b);
        assertOr(a.not(), b);
    }
    //Multi-Literal constructs
    public Lit and(Collection<Lit> args) {
        validate(args);
        return this.toLit(MonosatJNI.Ands(this.solverPtr, this.getLitBuffer(args), args.size()));
    }

    public Lit or(Collection<Lit> args) {
        validate(args);
        return this.toLit(MonosatJNI.Ors(this.solverPtr, this.getLitBuffer(args), args.size()));
    }

    public Lit nand(Collection<Lit> args) {
        validate(args);
        return this.toLit(MonosatJNI.Nands(this.solverPtr, this.getLitBuffer(args), args.size()));
    }

    public Lit nor(Collection<Lit> args) {
        validate(args);
        return this.toLit(MonosatJNI.Nors(this.solverPtr, this.getLitBuffer(args), args.size()));
    }

    public Lit xor(Collection<Lit> args) {
        validate(args);
        return this.toLit(MonosatJNI.Xors(this.solverPtr, this.getLitBuffer(args), args.size()));
    }

    public Lit xnor(Collection<Lit> args) {
        validate(args);
        return this.toLit(MonosatJNI.Xnors(this.solverPtr, this.getLitBuffer(args), args.size()));
    }
    public Lit implies(Lit a, Lit b){
        validate(a,b);
        return or(a.not(), b);
    }

    //assertion forms
    public void assertAnd(Collection<Lit> args) {
        validate(args);
        MonosatJNI.AssertAnds(this.solverPtr, this.getLitBuffer(args), args.size());
    }

    public void assertOr(Collection<Lit> args) {
        validate(args);
        MonosatJNI.AssertOrs(this.solverPtr, this.getLitBuffer(args), args.size());
    }

    public void assertNand(Collection<Lit> args) {
        validate(args);
        MonosatJNI.AssertNands(this.solverPtr, this.getLitBuffer(args), args.size());
    }

    public void assertNor(Collection<Lit> args) {
        validate(args);
        MonosatJNI.AssertNors(this.solverPtr, this.getLitBuffer(args), args.size());
    }

    public void assertXor(Collection<Lit> args) {
        validate(args);
        MonosatJNI.AssertXors(this.solverPtr, this.getLitBuffer(args), args.size());
    }

    public void assertXnor(Collection<Lit> args) {
        validate(args);
        MonosatJNI.AssertXnors(this.solverPtr, this.getLitBuffer(args), args.size());
    }

    //Bitvector constructs
    public BitVector ite(Lit condition, BitVector then, BitVector els) {
        validate(condition);
        validate(then,els);
        assert (then.width() == els.width());
        BitVector result = new BitVector(this, then.width());
        MonosatJNI.bv_ite(this.solverPtr, this.bvPtr, condition.l, then.id, els.id, result.id);
        return result;
    }

    public BitVector and(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_and(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public BitVector or(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_or(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public BitVector not(BitVector a) {
        validate(a);
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_not(solverPtr, bvPtr, a.id, result.id);
        return result;
    }

    public BitVector nand(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_nand(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public BitVector nor(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_nor(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public BitVector xor(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_xor(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public BitVector xnor(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_xnor(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public BitVector add(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_addition(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public BitVector subtract(BitVector a, BitVector b) {
        validate(a,b);
        assert(a.width()==b.width());
        BitVector result = new BitVector(this, a.width());
        MonosatJNI.bv_subtraction(solverPtr, bvPtr, a.id, b.id, result.id);
        return result;
    }

    public void assertEqual(BitVector a, BitVector b) {
        validate(a,b);
        assertTrue(a.geq(b));
        assertTrue(a.leq(b));
    }

    public void assertEqual(BitVector a, long constant) {
        validate(a);
        BitVector b = this.bv(a.width(), constant);
        assertTrue(a.geq(b));
        assertTrue(a.leq(b));
    }

    public void assertEqual(long constant, BitVector a) {
        validate(a);
        BitVector b = this.bv(a.width(), constant);
        assertTrue(a.geq(b));
        assertTrue(a.leq(b));
    }

    public BitVector min(Collection<BitVector> args) {
        validateBV(args);
        int w = args.iterator().next().width();
        BitVector result = new BitVector(this, w);
        MonosatJNI.bv_min(this.solverPtr, this.bvPtr, this.getBVBuffer(args, 0), args.size(), result.id);
        return result;
    }

    public BitVector min(BitVector a, BitVector b) {
        validate(a,b);
        ArrayList<BitVector> pair = new ArrayList<>();
        pair.add(a);
        pair.add(b);
        return min(pair);
    }

    public BitVector max(Collection<BitVector> args) {
        validateBV(args);
        int w = args.iterator().next().width();
        BitVector result = new BitVector(this, w);
        MonosatJNI.bv_min(this.solverPtr, this.bvPtr, this.getBVBuffer(args, 0), args.size(), result.id);
        return result;
    }

    public BitVector max(BitVector a, BitVector b) {
        validate(a,b);
        ArrayList<BitVector> pair = new ArrayList<>();
        pair.add(a);
        pair.add(b);
        return max(pair);
    }
}
