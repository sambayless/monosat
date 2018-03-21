import monosat.MonosatJNI;

import java.util.ArrayList;

public class BitVector {
    protected int id;
    int width;
    Solver solver;
    ArrayList<Lit> bits;

    public BitVector(Solver solver, ArrayList<Lit> bits){
        this.solver = solver;
        id = MonosatJNI.newBitvector(solver.solverPtr,solver.bvPtr,solver.getLitBuffer(bits),bits.size());
        width = bits.size();
        for(Lit l:bits){
            bits.add(l);
        }
    }

    public BitVector(Solver solver,int width, long constant){
        this.solver = solver;
        id = MonosatJNI.newBitvector_const(solver.solverPtr,solver.bvPtr,width,constant);
        this.width = width;
        for(int i = 0;i<width;i++){
            if ((constant & (1<<i))==1){
                bits.add(solver.True());
            }else{
                bits.add(solver.False());
            }
        }
    }

    /**
     * Creates a new BitVector if the specified width.
     * If introduceLiterals is true, then introduces width literals,
     * else creates a bitvector that has no literals associated with it.
     * @param solver
     * @param width
     */
    public BitVector(Solver solver, int width, boolean introduceLiterals){
        this.solver = solver;
        this.width = width;
        if (!introduceLiterals) {
            id = MonosatJNI.newBitvector_anon(solver.solverPtr, solver.bvPtr, width);
        }else{
            for(int i = 0;i<width;i++){
                bits.add(solver.newLit());
            }
            id = MonosatJNI.newBitvector(solver.solverPtr,solver.bvPtr,solver.getLitBuffer(bits),bits.size());
        }
    }

    public BitVector(Solver solver, int width){
        this(solver, width,true);
    }
    int width(){
        return width;
    }

    int newBVComparison_const_lt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
    int newBVComparison_bv_lt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
    int newBVComparison_const_leq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
    int newBVComparison_bv_leq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
    int newBVComparison_const_gt(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
    int newBVComparison_bv_gt(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);
    int newBVComparison_const_geq(SolverPtr S, BVTheoryPtr bv, int bvID, Weight weight);
    int newBVComparison_bv_geq(SolverPtr S, BVTheoryPtr bv, int bvID, int compareID);


    void bv_concat( SolverPtr S, BVTheoryPtr bv,int aID, int bID, int resultID);
    void bv_slice( SolverPtr S, BVTheoryPtr bv,int aID, int lower, int upper, int resultID);
    void bv_not( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvResultID);
    void bv_and( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
    void bv_nand( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
    void bv_or( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
    void bv_nor( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
    void bv_xor( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);
    void bv_xnor( SolverPtr S, BVTheoryPtr bv,int bvaID, int bvbID, int bvResultID);

    void bv_ite( SolverPtr S, BVTheoryPtr bv, int condition_lit,int bvThenID, int bvElseID, int bvResultID);

    void bv_addition( SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID);
    void bv_subtraction( SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID);
    void bv_multiply(SolverPtr S, BVTheoryPtr bv, int bvID1, int bvID2, int resultID);
    void bv_divide(SolverPtr S, BVTheoryPtr bv, int bvID1,  int bvID2, int resultID);
    void bv_min(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args,int resultID);
    void bv_max(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args,int resultID);
    void bv_popcount(SolverPtr S, BVTheoryPtr bv,  int* args,int n_args, int resultID);
    void bv_unary(SolverPtr S, BVTheoryPtr bv, int * args, int n_args, int resultID);


}
