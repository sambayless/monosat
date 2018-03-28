package monosat;
import monosat.MonosatJNI;
import monosat.Lit;
import monosat.Solver;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
public class BitVector {
    protected int id;
    int width;
    Solver solver;
    ArrayList<Lit> bits = new ArrayList<Lit>();

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
        assert(constant>=0);
        assert(constant<(1L<<width));
        id = MonosatJNI.newBitvector_const(solver.solverPtr,solver.bvPtr,width,constant);
        this.width = width;
        for(int i = 0;i<width;i++){
            if ((constant & (1<<i))==1){
                bits.add(solver.getTrue());
            }else{
                bits.add(solver.getFalse());
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
        if(width<0){
            throw new IllegalArgumentException("BitVector width must be >=0");
        }
        this.width = width;
        if (!introduceLiterals) {
            id = MonosatJNI.newBitvector_anon(solver.solverPtr, solver.bvPtr, width);
        }else{
            for(int i = 0;i<width;i++){
                bits.add(solver.newLit());
            }
            id = MonosatJNI.newBitvector(solver.solverPtr,solver.bvPtr,solver.getVarBuffer(bits,0),bits.size());
        }
    }

    public BitVector(Solver solver, int width){
        this(solver, width,true);
    }
    public List<Lit> getBits(){
        return Collections.unmodifiableList(bits);
    }

    public Lit get(int bit){
        return bits.get(bit);
    }

    int width(){
        return width;
    }
    int size(){
        return width();
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit gt(BitVector compareTo){
        int l = MonosatJNI.newBVComparison_bv_gt(solver.solverPtr,solver.bvPtr,this.id, compareTo.id);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit geq(BitVector compareTo){
        int l = MonosatJNI.newBVComparison_bv_geq(solver.solverPtr,solver.bvPtr,this.id, compareTo.id);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit lt(BitVector compareTo){
        int l = MonosatJNI.newBVComparison_bv_lt(solver.solverPtr,solver.bvPtr,this.id, compareTo.id);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit leq(BitVector compareTo){
        int l = MonosatJNI.newBVComparison_bv_leq(solver.solverPtr,solver.bvPtr,this.id, compareTo.id);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit neq(BitVector compareTo){
        int l = MonosatJNI.newBVComparison_bv_neq(solver.solverPtr,solver.bvPtr,this.id, compareTo.id);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit eq(BitVector compareTo){
        int l = MonosatJNI.newBVComparison_bv_eq(solver.solverPtr,solver.bvPtr,this.id, compareTo.id);
        return solver.toLit(l);
    }

    //Constant comparisons

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit gt(int compareTo){
        int l = MonosatJNI.newBVComparison_const_gt(solver.solverPtr,solver.bvPtr,this.id, compareTo);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit geq(int compareTo){
        int l = MonosatJNI.newBVComparison_const_geq(solver.solverPtr,solver.bvPtr,this.id, compareTo);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit lt(int compareTo){
        int l = MonosatJNI.newBVComparison_const_lt(solver.solverPtr,solver.bvPtr,this.id, compareTo);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit leq(int compareTo){
        int l = MonosatJNI.newBVComparison_const_leq(solver.solverPtr,solver.bvPtr,this.id, compareTo);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit neq(int compareTo){
        int l = MonosatJNI.newBVComparison_const_neq(solver.solverPtr,solver.bvPtr,this.id, compareTo);
        return solver.toLit(l);
    }
    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     * @param compareTo
     * @return
     */
    Lit eq(int compareTo){
        int l = MonosatJNI.newBVComparison_const_eq(solver.solverPtr,solver.bvPtr,this.id, compareTo);
        return solver.toLit(l);
    }

    /**
     * Creates a new bitvector consisting of the bits [this[0],..,this[size-1],append[0],..,append[append.size()-1]]
     * Does not introduce any new literals.
     * @param append BitVector to append.
     * @return
     */
    BitVector append(BitVector append){
        int w = width()+append.width();
        BitVector result = new BitVector(solver,w);
        MonosatJNI.bv_concat(solver.solverPtr,solver.bvPtr,this.id,append.id,result.id);
        return result;
    }

    /**
     * Create a new bitvector consisting of the bits [this[lower],..,this[upper-1]]
     * @param lower
     * @param upper
     * @return
     */
    BitVector slice(int lower, int upper){
        int w = upper-lower;
        assert(w>=0);
        BitVector result = new BitVector(solver,w);
        MonosatJNI.bv_slice(solver.solverPtr,solver.bvPtr,this.id,lower,upper-1,result.id);
        return result;
    }
    /**
     * Return a new bitvector consisting of the bitwise NOT of this bitvector.
     */
    BitVector not(){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_not(solver.solverPtr,solver.bvPtr,this.id,result.id);
        return result;
    }
    /**
     * Return a new bitvector consisting of the bitwise AND of this bitvector and other.
     */
    BitVector and(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_and(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }
    /**
     * Return a new bitvector consisting of the bitwise NAND of this bitvector and other.
     */
    BitVector nand(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_nand(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }
    /**
     * Return a new bitvector consisting of the bitwise OR of this bitvector and other.
     */
    BitVector or(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_or(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }
    /**
     * Return a new bitvector consisting of the bitwise NOR of this bitvector and other.
     */
    BitVector nor(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_nor(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }
    /**
     * Return a new bitvector consisting of the bitwise XOR of this bitvector and other.
     */
    BitVector xor(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_xor(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }
    /**
     * Return a new bitvector consisting of the bitwise XNOR of this bitvector and other.
     */
    BitVector xnor(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_xnor(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }

    /**
     * Returns a Bitvector that represents the non-wrapping two's complement addition
     * of this and other. To prevent wrapping, the solver will enforce that a+b<(1<<width()).
     * @param other
     * @return
     */
    BitVector add(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_addition(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }


    /**
     * Returns a Bitvector that represents the non-wrapping two's complement subtraction
     * of this and other. To prevent wrapping, the solver will enforce that a-b>=0.
     * @param other
     * @return
     */
    BitVector subtract(BitVector other){
        BitVector result = new BitVector(solver,width());
        MonosatJNI.bv_subtraction(solver.solverPtr,solver.bvPtr,this.id,other.id,result.id);
        return result;
    }




}
