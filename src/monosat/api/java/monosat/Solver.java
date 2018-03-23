import monosat.MonosatJNI;

import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Collection;

public class Solver {

    protected long solverPtr=0; //handle to the underlying monsoat solver instance.
    protected long bvPtr=0;
    int buffer_size0=1024;
    int buffer_size1=1024;
    int buffer_size2=1024;
    private IntBuffer ints0;
    private IntBuffer ints1;
    private IntBuffer ints2;

    Lit true_lit;
    //Instantiate a new solver
    public Solver(){
        solverPtr = MonosatJNI.newSolver();
        assert(solverPtr!=0);
        true_lit = new Lit(MonosatJNI.true_lit(solverPtr));
        initBV();
        initBuffers();
    }
    public Solver(String args){
        solverPtr = MonosatJNI.newSolver(args);
        assert(solverPtr!=0);
        initBV();
        initBuffers();
    }
    public Solver(ArrayList<String> args){
        String arg = "";
        for(String s:args){
            arg+= s + " ";
        }
        solverPtr = MonosatJNI.newSolver(arg);
        assert(solverPtr!=0);
        initBV();
        initBuffers();
    }
    private void  initBuffers(){
        initBuffer(0);
        initBuffer(1);
        initBuffer(2);
    }
    private void  initBuffer(int bufferN){
        if (bufferN==0){
            ints0 = ByteBuffer.allocateDirect(buffer_size0).asIntBuffer();
        }else if (bufferN==1){
            ints1 = ByteBuffer.allocateDirect(buffer_size1).asIntBuffer();
        }else if (bufferN==2) {
            ints2 = ByteBuffer.allocateDirect(buffer_size2).asIntBuffer();
        }
    }

    private IntBuffer  getBuffer(int bufferN, int minsize){

        if (bufferN==0){
            if (minsize>=buffer_size0){
                buffer_size0 = minsize*2;
                initBuffer(bufferN);
            }
            return ints0;
        }else if (bufferN==1){
            if (minsize>=buffer_size1){
                buffer_size1 = minsize*2;
                initBuffer(bufferN);
            }
            return ints1;
        }else if (bufferN==2) {
            if (minsize>=buffer_size2){
                buffer_size2 = minsize*2;
                initBuffer(bufferN);
            }
            return ints2;
        }
        throw new RuntimeException("BufferN must be between 0 and 3");
    }

    /**
     * Instantiate a bitvector theory solver
     */
    private void initBV(){
        assert(bvPtr==0);
        bvPtr = MonosatJNI.initBVTheory(solverPtr);
        assert(bvPtr!=0);
    }

    public void setOutputFile(String file){
        MonosatJNI.setOutputFile(solverPtr,file);
    }
    public Lit newLit(){
        return newLit(true);
    }
    public Lit newLit(boolean decisionVar){
        int var = MonosatJNI.newVar(solverPtr);
        assert(var>=0);
        Lit l = new Lit(var);
        while(lits.size()<=var){
            lits.add(null);
        }
        assert(!l.sign());
        lits.set(l.toVar(),l);
        assert(toLit(var*2)==l);
        return l;
    }

    void releaseLiteral(Lit l){
        l.validate();
        int x = l.l;
        assert(x>=0);
        MonosatJNI.releaseLiteral(solverPtr,x);
    }

    void setDecisionLiteral(Lit l, boolean decidable){
        MonosatJNI.setDecisionVar(solverPtr,l.toVar(), decidable);
    }
    boolean isDecisionLiteral(Lit l){
        return MonosatJNI.isDecisionVar(solverPtr,l.toVar());
    }
    public void setDecisionPriority(Lit l, int priority){
        MonosatJNI.setDecisionPriority(solverPtr, l.toVar(),priority);
    }
    public int getDecisionPriority(Lit l){
        return MonosatJNI.getDecisionPriority(solverPtr, l.toVar());
    }
    public void setDecisionPolarity(Lit l, boolean b){
        MonosatJNI.setDecisionPolarity(solverPtr, l.toVar(),b);
    }
    public boolean getDecisionPolarity(Lit l){
        return MonosatJNI.getDecisionPolarity(solverPtr, l.toVar());
    }
    public Lit True(){
        return true_lit;
    }
    public Lit False() {
        return true_lit.negate();
    }
    public void disallowSimplification(Lit l){
        MonosatJNI.disallowLiteralSimplification(solverPtr, l.toVar());
    }
    public void disablePreprocessing(){
        MonosatJNI.disablePreprocessing(solverPtr);
    }
    public int nVars(){
        return MonosatJNI.nVars(solverPtr);
    }
    public int nClauses(){
        return MonosatJNI.nClauses(solverPtr);
    }
    public int nBitvectors(){
        return MonosatJNI.nBitvectors(solverPtr,bvPtr);
    }

    IntBuffer getBVBuffer(Collection<BitVector> clause,int bufferN){
        assert(bufferN<3);
        assert(bufferN>=0);
        IntBuffer buffer = getBuffer(bufferN, clause.size());
        int index =0;
        for(BitVector bv:clause){

            buffer.put(index,bv.id);
            index++;
        }
        return buffer;
    }

    IntBuffer getLitBuffer(Collection<Lit> clause){
        return getLitBuffer(clause,0);
    }
    IntBuffer getLitBuffer(Lit[] clause,int bufferN){
        assert(bufferN<3);
        assert(bufferN>=0);
        IntBuffer buffer = getBuffer(bufferN, clause.length);
        int index =0;
        for(Lit l:clause){
            l.validate();
            buffer.put(index,l.l);
            index++;
        }
        return buffer;
    }
    IntBuffer getLitBuffer(Collection<Lit> clause,int bufferN){
        assert(bufferN<3);
        assert(bufferN>=0);
        IntBuffer buffer = getBuffer(bufferN, clause.size());
        int index =0;
        for(Lit l:clause){
            l.validate();
            buffer.put(index,l.l);
            index++;
        }
        return buffer;
    }
    IntBuffer getIntBuffer(Collection<Integer> ints,int bufferN){
        assert(bufferN<3);
        assert(bufferN>=0);
        IntBuffer buffer = getBuffer(bufferN, ints.size());
        int index =0;
        for(Integer i:ints){
            buffer.put(index,i);
            index++;
        }
        return buffer;
    }
    /**
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     * @param clause The clause to add to the solver
     * @return
     */
    public boolean addClause(ArrayList<Lit> clause){

        boolean status = MonosatJNI.addClause(solverPtr,getLitBuffer(clause,1),clause.size());
        return status;
    }
    /**
     * Returns false if the formula is trivially unsatisfiable after adding this clause (or if it was already trivially unsatisfiable),
     * else returns true.
     * @param l The literal to add as a unit clause to the solver (forcing l to be true)
     * @return
     */
    public boolean addClause(Lit l){
        l.validate();
        boolean status = MonosatJNI.addUnitClause(solverPtr,l.l);
        return status;
    }
    public boolean addClause(Lit l1, Lit l2){
        l1.validate();l2.validate();
        boolean status = MonosatJNI.addBinaryClause(solverPtr,l1.l, l2.l);
        return status;
    }
    public boolean addClause(Lit l1, Lit l2, Lit l3){
        l1.validate();l2.validate();l3.validate();
        boolean status = MonosatJNI.addTertiaryClause(solverPtr,l1.l, l2.l, l3.l);
        return status;
    }

    //Solver API

    //basic solver functions
    public boolean solve(){
        boolean r = MonosatJNI.solve(solverPtr);
        return r;
    }
    public boolean solve(Lit... assumptions){
        return MonosatJNI.solveAssumptions(solverPtr, getLitBuffer(assumptions,0),assumptions.length);
    }

    public boolean solve(Collection<Lit> assumptions){
        return MonosatJNI.solveAssumptions(solverPtr, getLitBuffer(assumptions),assumptions.size());
    }



    //Sets the (approximate) time limit in seconds before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable time limit.
    public void setTimeLimit(int seconds){
        MonosatJNI.setTimeLimit(solverPtr,seconds);
    }
    //Sets the (approximate) memory limit in megabytes before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable memory limit.
    public void setMemoryLimit(int mb){
        MonosatJNI.setMemoryLimit(solverPtr, mb);
    }
    //Sets the maximum number of (additional) conflicts allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable conflict limit.
    public void setConflictLimit(int num_conflicts){
        MonosatJNI.setConflictLimit(solverPtr,num_conflicts);
    }
    //Sets the maximum number of (additional) propagations allowed in the solver before returning l_Undef from solveLimited; ignored by solve(). Set to <0 to disable propagation limit.
    public void setPropagationLimit(int num_propagations){
        MonosatJNI.setPropagationLimit(solverPtr,num_propagations);
    }

    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    //Consider using Optional<Boolean> instead.
    public LBool solveLimited(){
        int result = MonosatJNI.solveLimited(solverPtr);
        assert(result>=0);
        assert(result<=2);
        return LBool.values()[result];
    }
    //Returns 0 for satisfiable, 1 for proved unsatisfiable, 2 for failed to find a solution (within any resource limits that have been set)
    public LBool solveAssumptionsLimited(ArrayList<Lit> assumptions){
        int result = MonosatJNI.solveAssumptionsLimited(solverPtr,getLitBuffer(assumptions),assumptions.size());
        assert(result>=0);
        assert(result<=2);
        return LBool.values()[result];
    }

    public boolean lastSolutionWasOptimal(){
        return MonosatJNI.lastSolutionWasOptimal(solverPtr);
    }

    /**
     * If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
     */
    public ArrayList<Lit>  getConflictClause(ArrayList<Lit> store){
        store.clear();

        //If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a subset of the assumptions which are sufficient to cause the instance to be UNSAT).
        //Fills the given pointer with the first max_store_size literals of the conflict clause, and returns the number of literals in the conflict clause. Set store_clause to null and max_store_size to 0 to find the size of the conflict clause
        //Returns -1 if the solver has no conflict clause from the most recent solve() call (because that call was not UNSAT)
        int conflict_size = MonosatJNI.getConflictClause(solverPtr,null,0);
        IntBuffer buf = getBuffer(0,conflict_size);
        int sz = MonosatJNI.getConflictClause(solverPtr, buf, conflict_size);
        assert(sz == conflict_size);
        for(int i = 0;i<conflict_size;i++){
            int l = buf.get(i);
            Lit lit = toLit(l);
            store.add(lit);
        }
        return store;
    }
    public ArrayList<Lit>  getConflictClause(){
        ArrayList<Lit> store = new ArrayList<Lit>();
        getConflictClause(store);
        return store;
    }

    //Backtrack the solver to level 0
    void backtrack(){
        MonosatJNI.backtrack(solverPtr);
    }


    //Holds positive versions of all literals, so that we don't need to create multiple literal objects for the same literal
    //As literals have no reference to the solver and are really just a thin wrapper around an integer,
    //this list can be safely shared across all solvers.
    private static ArrayList<Lit> lits;

    protected Lit toLit(int literal){
        assert(literal>=0);
        int var = literal/2;
        assert(var<nVars());//the variable must have already been declared in the sat solver before this call
        while(var>=lits.size()){
            lits.add(null);
        }
        if(lits.get(var)==null){
            lits.set(var,new Lit(var));
        }
        Lit l = lits.get(var);
        if ((literal &1) ==1){
            return l.negate();
        }else{
            return l;
        }
    }

    //Optimization API

    public void clearOptimizationObjectives(){
        MonosatJNI.clearOptimizationObjectives(solverPtr);
    }

    public void maximizeBV(BitVector bv){
        MonosatJNI.maximizeBV(solverPtr,bvPtr,bv.id);
    }


    public void minimizeBV(BitVector bv){
        MonosatJNI.minimizeBV(solverPtr,bvPtr,bv.id);
    }
    public void maximizeLits(ArrayList<Lit> literals){
        MonosatJNI.maximizeLits(solverPtr,getLitBuffer(literals),literals.size());
    }
    public void minimizeLits(ArrayList<Lit> literals){
        MonosatJNI.minimizeLits(solverPtr,getLitBuffer(literals),literals.size());
    }

    public void maximizeWeightedLits(ArrayList<Lit> literals, ArrayList<Integer> weights){
        assert(literals.size()==weights.size());
        MonosatJNI.maximizeWeightedLits(solverPtr,getLitBuffer(literals), getIntBuffer(weights,1), literals.size());
    }
    public void minimizeWeightedLits(ArrayList<Lit> literals, ArrayList<Integer> weights){
        assert(literals.size()==weights.size());
        MonosatJNI.minimizeWeightedLits(solverPtr,getLitBuffer(literals), getIntBuffer(weights,1), literals.size());
    }

    void AssertAtMostOne(ArrayList<Lit> clause){
        MonosatJNI.at_most_one(solverPtr,getLitBuffer(clause),clause.size());
    }
    enum Comparison{LT,LEQ,EQ,GEQ,GT}

    void AssertPB(ArrayList<Lit> clause, int compareTo,Comparison c){
        AssertPB(clause,null,compareTo,c);
    }
    void AssertPB(ArrayList<Lit> clause,ArrayList<Integer> weights, int compareTo,Comparison c){
        IntBuffer wt_buffer = getBuffer(1,clause.size());
        int n_wts =0;
        if(weights!=null) {
            n_wts = Math.min(clause.size(),weights.size());
            for (int i = 0; i < clause.size() && i < weights.size(); i++) {
                wt_buffer.put(i,weights.get(i));
            }
        }
        for(int i = n_wts;i<clause.size();i++){
            wt_buffer.put(i,1); //default weight is 1
        }
        switch(c){
            case LT:
                MonosatJNI.assertPB_lt(solverPtr,compareTo,clause.size(),getLitBuffer(clause),wt_buffer);
            case LEQ:
                MonosatJNI.assertPB_leq(solverPtr,compareTo,clause.size(),getLitBuffer(clause),wt_buffer);
            case EQ:
                MonosatJNI.assertPB_eq(solverPtr,compareTo,clause.size(),getLitBuffer(clause),wt_buffer);
            case GEQ:
                MonosatJNI.assertPB_geq(solverPtr,compareTo,clause.size(),getLitBuffer(clause),wt_buffer);
            case GT:
                MonosatJNI.assertPB_gt(solverPtr,compareTo,clause.size(),getLitBuffer(clause),wt_buffer);
        }
    }

    //Convert any pb constraints in the solver into cnf (will be called automatically before solve())
    public void flushPB(){
        MonosatJNI.flushPB(solverPtr);
    }

    private ArrayList<ArrayList<BitVector>> cached_bvs;
    private static long MAX_CACHE_CONST = 255;
    /**
     * Caches constant bitvectors of value < MAX_CACHE_CONST
     * @param width
     * @param constant
     * @return
     */
    public BitVector bv(int width, long constant){
        if(constant>=0 && constant<=MAX_CACHE_CONST){
            while(cached_bvs.size()<=width){
                cached_bvs.add(new ArrayList<>());
            }

            int small_const = (int) constant;
            while(cached_bvs.get(width).size()<=constant){
                cached_bvs.get(width).add(null);
            }
            if(cached_bvs.get(width).get(small_const)==null){
                cached_bvs.get(width).set(small_const,new BitVector(this,width,constant));
            }
            assert(cached_bvs.get(width).get(small_const)!=null);
            return cached_bvs.get(width).get(small_const);
        }
        return new BitVector(this,width,constant);
    }



}
