import monosat.MonosatJNI;

import java.util.ArrayList;
import java.util.Collection;

public class Logic {
    private static ThreadLocal<Solver> _solver;
    public static Solver getSolver(){
        return _solver.get();
    }

 
    //Literal level constructs

    public static Lit ite(Lit condition, Collection<Lit> then,Collection<Lit> els){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Ite(solver.solverPtr,condition.l,then.l,els.l));
    }

    public static Lit and(Lit a, Lit b){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.And(solver.solverPtr,a.l,b.l));
    }
    public static Lit or(Lit a, Lit b){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Or(solver.solverPtr,a.l,b.l));
    }
    public static Lit not(Lit a){
        return a.negate();
    }
    public static Lit nand(Lit a, Lit b){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Nand(solver.solverPtr,a.l,b.l));
    }
    public static Lit nor(Lit a, Lit b){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Nor(solver.solverPtr,a.l,b.l));
    }
    public static Lit xor(Lit a, Lit b){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Xor(solver.solverPtr,a.l,b.l));
    }
    public static Lit xnor(Lit a, Lit b){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Xnor(solver.solverPtr,a.l,b.l));
    }

    //Assertion forms. Note, these are capitalized to avoid name-classing with the builtin 'Assert' method
    public static void Assert(Lit a){
        MonosatJNI.Assert(getSolver().solverPtr,a.l);
    }
    public static void AssertNot(Lit a){
        MonosatJNI.Assert(getSolver().solverPtr,a.negate().l);
    }

    public static void AssertAnd(Lit a, Lit b){
        Solver solver = getSolver();
        MonosatJNI.AssertAnd(solver.solverPtr,a.l,b.l);
    }
    public static void AssertOr(Lit a, Lit b){
        Solver solver = getSolver();
       MonosatJNI.AssertOr(solver.solverPtr,a.l,b.l);
    }

    public static void AssertNand(Lit a, Lit b){
        Solver solver = getSolver();
        MonosatJNI.Nand(solver.solverPtr,a.l,b.l);
    }
    public static void AssertNor(Lit a, Lit b){
        Solver solver = getSolver();
        MonosatJNI.AssertNor(solver.solverPtr,a.l,b.l);
    }
    public static void AssertXor(Lit a, Lit b){
        Solver solver = getSolver();
        MonosatJNI.AssertXor(solver.solverPtr,a.l,b.l);
    }
    public static void AssertXnor(Lit a, Lit b){
        Solver solver = getSolver();
        MonosatJNI.AssertXnor(solver.solverPtr,a.l,b.l);
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
        return solver.toLit(MonosatJNI.And(solver.solverPtr,a.l,b.l));
    }
    public static BitVector or(BitVector a, BitVector b){
        Solver solver = getSolver();
        return solver.toLit(MonosatJNI.Or(solver.solverPtr,a.l,b.l));
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
        AssertNot(a.gt(b));
        AssertNot(a.lt(b));
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
