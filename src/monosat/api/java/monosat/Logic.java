package monosat;
import monosat.MonosatJNI;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

public final class Logic {
    private static final ThreadLocal<Solver> _solver = new ThreadLocal();

    private Logic(){

    }

    public static Solver newSolver(){
        setSolver(new Solver());
        return getSolver();
    }

    public static Solver newSolver(String args){
        setSolver(new Solver(args));
        return getSolver();
    }
    public static void setSolver(Solver s){
        _solver.set(s);
    }

    public static Solver getSolver(){
        if (_solver.get() == null){
            newSolver();
        }
        return _solver.get();
    }

    public static Lit newLit(){
        return newLit(true);
    }
    public static Lit newLit(boolean decidable){
        return getSolver().newLit(decidable);
    }

    public static Lit True(){
        return getSolver().True();
    }
    public static Lit False(){
        return getSolver().False();
    }
    public static boolean solve(){
        return getSolver().solve();
    }

    public static boolean solve(Lit...assumptions){
        return getSolver().solve(assumptions);
    }

    public static boolean solve(Collection<Lit> assumptions){
        return getSolver().solve(assumptions);
    }
    //Literal level constructs
    public static Lit ite(Lit condition, Lit then,Lit els){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Ite(solver.solverPtr,condition.l,then.l,els.l));
    }

    public static Lit and(Lit ... args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            return args[0];
        }else if (args.length==2){
            Solver solver = getSolver();
            return solver.toLit(MonosatJNI.And(solver.solverPtr,args[0].l,args[1].l));
        }
        return and(Arrays.asList(args));
    }
    public static Lit or(Lit... args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            return args[0];
        }else if (args.length==2){
            Solver solver = getSolver();
            return solver.toLit(MonosatJNI.Or(solver.solverPtr,args[0].l,args[1].l));
        }
        return or(Arrays.asList(args));
    }
    public static Lit not(Lit a){
        return a.negate();
    }
    public static Lit nand(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            return args[0].negate();
        }else if (args.length==2){
            Solver solver = getSolver();
            return solver.toLit(MonosatJNI.Nand(solver.solverPtr,args[0].l,args[1].l));
        }
        return nand(Arrays.asList(args));
    }
    public static Lit nor(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            return args[0];
        }else if (args.length==2){
            Solver solver = getSolver();
            return solver.toLit(MonosatJNI.Nor(solver.solverPtr,args[0].l,args[1].l));
        }
        return nor(Arrays.asList(args));
    }
    public static Lit xor(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            return args[0].negate();
        }else if (args.length==2){
            Solver solver = getSolver();
            return solver.toLit(MonosatJNI.Xor(solver.solverPtr,args[0].l,args[1].l));
        }
        return xor(Arrays.asList(args));
    }
    public static Lit xnor(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            return args[0];
        }else if (args.length==2){
            Solver solver = getSolver();
            return solver.toLit(MonosatJNI.Xnor(solver.solverPtr,args[0].l,args[1].l));
        }
        return xnor(Arrays.asList(args));
    }

    //Assertion forms. Note, these are capitalized to avoid name-classing with the builtin 'Assert' method
    public static void Assert(Lit a){
        MonosatJNI.Assert(getSolver().solverPtr,a.l);
    }
    public static void AssertNot(Lit a){
        MonosatJNI.Assert(getSolver().solverPtr,a.negate().l);
    }

    public static void AssertAnd(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            Assert(args[0]);
        }else if (args.length==2){
            MonosatJNI.AssertAnd(getSolver().solverPtr,args[0].l,args[1].l);
        }
         AssertAnd(Arrays.asList(args));
    }
    public static void AssertOr(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            Assert(args[0]);
        }else if (args.length==2){
            MonosatJNI.AssertOr(getSolver().solverPtr,args[0].l,args[1].l);
        }
        AssertOr(Arrays.asList(args));
    }

    public static void AssertNand(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            Assert(args[0].negate());
        }else if (args.length==2){
            MonosatJNI.AssertNand(getSolver().solverPtr,args[0].l,args[1].l);
        }
        AssertNand(Arrays.asList(args));
    }
    public static void AssertNor(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            Assert(args[0].negate());
        }else if (args.length==2){
            MonosatJNI.AssertNor(getSolver().solverPtr,args[0].l,args[1].l);
        }
        AssertNor(Arrays.asList(args));
    }
    public static void AssertXor(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            Assert(args[0]);//If this the correct behaviour here?
        }else if (args.length==2){
            MonosatJNI.AssertXor(getSolver().solverPtr,args[0].l,args[1].l);
        }
        AssertXor(Arrays.asList(args));
    }
    public static void AssertXnor(Lit...args){
        if (args.length==0) {
            throw new IllegalArgumentException("Requires at least 1 argument");
        }else if (args.length==1){
            Assert(args[0]);//If this the correct behaviour here?
        }else if (args.length==2){
            MonosatJNI.AssertXnor(getSolver().solverPtr,args[0].l,args[1].l);
        }
        AssertXnor(Arrays.asList(args));
    }
    public static void AssertEqual(Lit a, Lit b){
        AssertXnor(a,b);
    }

    //Multi-Literal constructs
    public static Lit and(Collection<Lit> elements){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Ands(solver.solverPtr,solver.getLitBuffer(elements),elements.size()));
    }
    public static Lit or(Collection<Lit> elements){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Ors(solver.solverPtr,solver.getLitBuffer(elements),elements.size()));
    }
    public static Lit nand(Collection<Lit> elements){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Nands(solver.solverPtr,solver.getLitBuffer(elements),elements.size()));
    }
    public static Lit nor(Collection<Lit> elements){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Nors(solver.solverPtr,solver.getLitBuffer(elements),elements.size()));
    }
    public static Lit xor(Collection<Lit> elements){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Xors(solver.solverPtr,solver.getLitBuffer(elements),elements.size()));
    }
    public static Lit xnor(Collection<Lit> elements){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Xnors(solver.solverPtr,solver.getLitBuffer(elements),elements.size()));
    }

    //Assertion forms
    public static void Assert(Collection<Lit> elements){
        AssertOr(elements);
    }
    public static void AssertAnd(Collection<Lit> elements){
        Solver solver = getSolver();
        MonosatJNI.AssertAnds(solver.solverPtr,solver.getLitBuffer(elements),elements.size());
    }
    public static void AssertOr(Collection<Lit> elements){
        Solver solver = getSolver();
        MonosatJNI.AssertOrs(solver.solverPtr,solver.getLitBuffer(elements),elements.size());
    }
    public static void AssertNand(Collection<Lit> elements){
        Solver solver = getSolver();
        MonosatJNI.AssertNands(solver.solverPtr,solver.getLitBuffer(elements),elements.size());
    }
    public static void AssertNor(Collection<Lit> elements){
        Solver solver = getSolver();
        MonosatJNI.AssertNors(solver.solverPtr,solver.getLitBuffer(elements),elements.size());
    }
    public static void AssertXor(Collection<Lit> elements){
        Solver solver = getSolver();
        MonosatJNI.AssertXors(solver.solverPtr,solver.getLitBuffer(elements),elements.size());
    }
    public static void AssertXnor(Collection<Lit> elements){
        Solver solver = getSolver();
        MonosatJNI.AssertXnors(solver.solverPtr,solver.getLitBuffer(elements),elements.size());
    }


    //Bitvector constructs
    public static BitVector ite(Lit condition, BitVector then, BitVector els){
        Solver solver = getSolver();
        assert(then.width()==els.width());
        BitVector result = new BitVector(solver, then.width());
        MonosatJNI.bv_ite(solver.solverPtr,solver.bvPtr,condition.l,then.id,els.id,result.id);
        return result;
    }
    public static BitVector and(BitVector a, BitVector b){
        Solver solver = getSolver();
        return a.and(b);
    }
    public static BitVector or(BitVector a, BitVector b){
        return a.or(b);
    }
    public static BitVector not(BitVector a){
        return a.not();
    }
    public static BitVector nand(BitVector a, BitVector b){
        return a.nand(b);
    }
    public static BitVector nor(BitVector a, BitVector b){
        return a.nor(b);
    }
    public static BitVector xor(BitVector a, BitVector b){
        return a.xor(b);
    }
    public static BitVector xnor(BitVector a, BitVector b){
        return a.xnor(b);
    }

    public static BitVector add(BitVector a, BitVector b){
        return a.add(b);
    }
    public static BitVector subtract(BitVector a, BitVector b){
        return a.subtract(b);
    }

    public static void AssertEqual(BitVector a, BitVector b){
        Assert(a.geq(b));
        Assert(a.leq(b));
    }

    public static void AssertEqual(BitVector a, long constant){
        BitVector b = getSolver().bv(a.width(),constant);
        Assert(a.geq(b));
        Assert(a.leq(b));
    }
    public static void AssertEqual(long constant,BitVector a){
        BitVector b = getSolver().bv(a.width(),constant);
        Assert(a.geq(b));
        Assert(a.leq(b));
    }
    public static BitVector min(Collection<BitVector> args){
        Solver solver = getSolver();
        assert(args.size()>=0);
        int w =args.iterator().next().width();
        BitVector result = new BitVector(solver, w);
        MonosatJNI.bv_min(solver.solverPtr,solver.bvPtr,solver.getBVBuffer(args,0),args.size(),result.id);
        return result;
    }
    public static BitVector min(BitVector a, BitVector b){
        ArrayList pair = new ArrayList();
        pair.add(a);
        pair.add(b);
        return min(pair);
    }
    public static BitVector max(Collection<BitVector> args){
        Solver solver = getSolver();
        assert(args.size()>=0);
        int w =args.iterator().next().width();
        BitVector result = new BitVector(solver, w);
        MonosatJNI.bv_min(solver.solverPtr,solver.bvPtr,solver.getBVBuffer(args,0),args.size(),result.id);
        return result;
    }
    public static BitVector max(BitVector a, BitVector b){
        ArrayList pair = new ArrayList();
        pair.add(a);
        pair.add(b);
        return max(pair);
    }

    /**
     * Create a new bitvector  of the specified width, with a constant value
     * @param width
     * @param constant
     * @return
     */
    public static BitVector bv(int width, long constant){
        return getSolver().bv(width,constant);
    }

    /**
     * Create a new bitvector of the specified width.
     * @param width
     * @return
     */
    public static BitVector bv(int width){
        return new BitVector(getSolver(),width);
    }
}
