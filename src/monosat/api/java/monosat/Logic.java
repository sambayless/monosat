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

import java.nio.IntBuffer;
import java.util.*;

/**
 * Logic provides static accessors for common logic functions.
 * <p>
 * The expected use of Logic is for all methods of Logic to be statically imported:
 * import static monosat.Logic.*
 * These static methods form a light-weight domain specific language.
 */
public final class Logic {
    //prevent instances of logic from being constructed
    private Logic() {}
    private static boolean allow_contradictions=false;

    private static Solver getSolver(Lit... args){
        for(Lit l:args){
            if(l==null || l==Lit.Error || l==Lit.Undef){
                throw new IllegalArgumentException("Invalid literal " + String.valueOf(l));
            }
            if (l.solver!=null){
                return l.solver;
            }
        }
        return null;
    }

    private static Solver getSolver(Collection<Lit> args){
        for(Lit l:args){
            if(l==null || l==Lit.Error || l==Lit.Undef){
                throw new IllegalArgumentException("Invalid literal " + String.valueOf(l));
            }
            if (l.solver!=null){
                return l.solver;
            }
        }
        return null;
    }

    /**
     * By default, certain trivial contradictions (such as AssertFalse(True), or AssertOr() with empty arguments)
     * will be caught, and trigger an exception. This is useful, as such contradictions are almost always unintentional.
     * Use this method to disable exceptions on trivial contradictions.
     */
    public static synchronized void allowContradictions(){
        allow_contradictions=true;
    }

    private static synchronized void contradiction(Lit a){
        //all existing solvers are now UNSAT
        for(Solver s: Solver.solvers.keySet()){
            s.assertFalse(Lit.True);
        }
        if(a==Lit.False){
            throw new RuntimeException("Constant False asserted to be True (which is almost always an error).");
        }else if(a==Lit.True){
            throw new RuntimeException("Constant True asserted to be False (which is almost always an error).");
        }
    }

    private static synchronized void contradiction(){
        //all existing solvers are now UNSAT
        for(Solver s: Solver.solvers.keySet()){
            s.assertFalse(Lit.True);
        }
        throw new RuntimeException("Statically UNSAT assertion (which is almost always an error).");
    }

    public static Lit not(Lit a) {
        return a.not();
    }

    public static Lit implies(Lit a, Lit b){
        return a.getSolver().implies(a, b);
    }

    public static Lit ite(Lit condition, Lit then, Lit els) {
        Solver solver = getSolver(condition,then,els);
        if(solver!=null) {
            return solver.ite(condition, then, els);
        }else{
            if(condition==Lit.True){
                return then;
            }else{
                return els;
            }
        }
    }

    public static Lit and(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.and(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    return Lit.False;
                }
            }
            return Lit.True;
        }
    }

    public static Lit and(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.and(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    return Lit.False;
                }
            }
            return Lit.True;
        }
    }

    public static Lit or(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.or(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    return Lit.True;
                }
            }
            return Lit.False;
        }
    }
    public static Lit or(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.or(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    return Lit.True;
                }
            }
            return Lit.False;
        }
    }


    public static Lit nand(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.nand(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    return Lit.True;
                }
            }
            return Lit.False;
        }
    }


    public static Lit nand(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.nand(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    return Lit.True;
                }
            }
            return Lit.False;
        }
    }

    public static Lit nor(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.nor(args);
        }else{
            for(Lit l:args){
                if(l!=Lit.False){
                    return Lit.False;
                }
            }
            return Lit.True;
        }
    }
    public static Lit nor(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.nor(args);
        }else{
            for(Lit l:args){
                if(l!=Lit.False){
                    return Lit.False;
                }
            }
            return Lit.True;
        }
    }
    public static Lit xor(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.xor(args);
        }else{
            //XOR is defined to be true if an odd number of its arguments are true, and false otherwise
           int trueCount = 0;
           for(Lit l:args){
               if(l==Lit.True){
                   trueCount++;
               }
           }
           if(trueCount%2==1){
               return Lit.True;
           }else{
               return Lit.False;
           }
        }
    }
    public static Lit xor(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.xor(args);
        }else{
            //XOR is defined to be true if an odd number of its arguments are true, and false otherwise
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            if(trueCount%2==1){
                return Lit.True;
            }else{
                return Lit.False;
            }
        }
    }
    public static Lit xnor(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.xnor(args);
        }else{
            //XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            if(trueCount%2==0){
                return Lit.True;
            }else{
                return Lit.False;
            }
        }
    }
    public static Lit xnor(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            return solver.xnor(args);
        }else{
            //XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            if(trueCount%2==0){
                return Lit.True;
            }else{
                return Lit.False;
            }
        }
    }
    public static Lit equal(Lit a, Lit b) {
        return xnor(a,b);
    }
    //Assertion forms.
    public static void assertTrue(Lit a) {
        Solver solver = getSolver(a);
        if(solver!=null){
            solver.assertTrue(a);
        }else{
            if(a==Lit.False){
                contradiction(a);
            }else{
                //do nothing
            }
        }
    }

    public static void assertFalse(Lit a) {
        Solver solver = getSolver(a);
        if(solver!=null){
            solver.assertFalse(a);
        }else{
            if(a==Lit.True){
                contradiction(a);
            }else{
                //do nothing
            }
        }
    }

    public static void assertAnd(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertAnd(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    contradiction(l);
                }
            }
            //do nothing
        }
    }
    public static void assertAnd(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertAnd(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    contradiction(l);
                }
            }
            //do nothing
        }
    }
    public static void assertOr(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertOr(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    return;
                }
            }
            contradiction();
        }
    }
    public static void assertOr(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertOr(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    return;
                }
            }
            contradiction();
        }
    }
    public static void assertNand(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertNand(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    return;
                }
            }
            contradiction();
        }
    }
    public static void assertNand(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertNand(args);
        }else{
            for(Lit l:args){
                if(l==Lit.False){
                    return;
                }
            }
            contradiction();
        }
    }
    public static void assertNor(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertNor(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    contradiction(l);
                }
            }
        }
    }
    public static void assertNor(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertNor(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    contradiction(l);
                }
            }
        }
    }
    public static void assertXor(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertXor(args);
        }else{
            //XOR is defined to be true if an odd number of its arguments are true, and false otherwise
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            if(trueCount%2==1){
                //do nothing
            }else{
                contradiction();
            }
        }
    }
    public static void assertXor(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertXor(args);
        }else{
            //XOR is defined to be true if an odd number of its arguments are true, and false otherwise
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            if(trueCount%2==1){
                //do nothing
            }else{
                contradiction();
            }
        }
    }
    public static void assertXnor(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertXnor(args);
        }else{
            //XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            if(trueCount%2==1){
                contradiction();
            }else{

            }
        }
    }
    public static void assertXnor(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertXnor(args);
        }else{
            //XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            if(trueCount%2==1){
                contradiction();
            }else{

            }
        }
    }
    public static void assertEqual(Lit a, Lit b) {
        Solver solver = getSolver(a,b);
        if(solver!=null) {
           solver.assertEqual(a, b);
        }else{
            if(a!=b){
                contradiction();
            }
        }
    }

    public static void assertImplies(Lit a, Lit b) {
        Solver solver = getSolver(a,b);
        if(solver!=null) {
            solver.assertImplies(a, b);
        }else{
            if(a==Lit.True && b==Lit.False){
                contradiction();
            }
        }
    }

    //Bitvector constructs
    public static BitVector ite(Lit condition, BitVector then, BitVector els) {
        return then.getSolver().ite(condition, then, els);
    }

    public static BitVector and(BitVector a, BitVector b) {
        return a.getSolver().and(a, b);
    }

    public static BitVector or(BitVector a, BitVector b) {
        return a.getSolver().or(a, b);
    }

    public static BitVector not(BitVector a) {
        return a.getSolver().not(a);
    }

    public static BitVector nand(BitVector a, BitVector b) {
        return a.getSolver().nand(a, b);
    }

    public static BitVector nor(BitVector a, BitVector b) {
        return a.getSolver().nor(a, b);
    }

    public static BitVector xor(BitVector a, BitVector b) {
        return a.getSolver().xor(a, b);
    }

    public static BitVector xnor(BitVector a, BitVector b) {
        return a.getSolver().xnor(a, b);
    }

    public static BitVector add(BitVector a, BitVector b) {
        return a.getSolver().add(a, b);
    }

    public static BitVector subtract(BitVector a, BitVector b) {
        return a.getSolver().subtract(a, b);
    }

    public static void assertEqual(BitVector a, BitVector b) {
        a.getSolver().assertEqual(a, b);
    }

    public static void assertEqual(BitVector a, long constant) {
        a.getSolver().assertEqual(a, constant);
    }

    public static void assertEqual(long constant, BitVector a) {
        a.getSolver().assertEqual(constant, a);
    }

    public static BitVector min(Collection<BitVector> args) {
        if(args.size()>0) {
            return args.iterator().next().getSolver().min(args);
        }else{
            throw new IllegalArgumentException("min requires at least one argument");
        }
    }

    public static BitVector min(BitVector a, BitVector b) {
        return a.getSolver().min(a, b);
    }

    public static BitVector max(Collection<BitVector> args) {
        if(args.size()>0) {
            return args.iterator().next().getSolver().max(args);
        }else{
            throw new IllegalArgumentException("max requires at least one argument");
        }
    }

    public static BitVector max(BitVector a, BitVector b) {
        return a.getSolver().max(a, b);
    }

    //Psuedo-Boolean constraints

    public static void maximizeBV(BitVector bv) {
        bv.getSolver().maximizeBV(bv);
    }

    public static void minimizeBV(BitVector bv) {
        bv.getSolver().minimizeBV(bv);
    }

    public static void maximizeLits(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.maximizeLits(args);
        }else{
            //do nothing
        }
    }

    public static void minimizeLits(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.minimizeLits(args);
        }else{
            //do nothing
        }
    }

    public static void maximizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        Solver solver = getSolver(literals);
        if(solver!=null){
            solver.maximizeWeightedLits(literals, weights);
        }else{
            //do nothing
        }
    }

    public static void minimizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        Solver solver = getSolver(literals);
        if(solver!=null){
            solver.minimizeWeightedLits(literals, weights);
        }else{
            //do nothing
        }
    }

    public static void assertAtMostOne(Lit... args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertAtMostOne(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    return;
                }
            }
            contradiction();
        }
    }

    public static void assertAtMostOne(Collection<Lit> args) {
        Solver solver = getSolver(args);
        if(solver!=null){
            solver.assertAtMostOne(args);
        }else{
            for(Lit l:args){
                if(l==Lit.True){
                    return;
                }
            }
            contradiction();
        }
    }
    private static boolean checkStaticPB(int trueCount, Comparison c, int compareTo){
        switch(c){
            case GT:
                return trueCount>compareTo;
            case GEQ:
                return trueCount>=compareTo;
            case EQ:
                return trueCount==compareTo;
            case LT:
                return trueCount<compareTo;
            case LEQ:
                return trueCount<=compareTo;
            case NEQ:
                return trueCount!=compareTo;
        }
        assert(false);
        return false;
    }

    public static void assertPB(Collection<Lit> args, Comparison c, int compareTo) {
        Solver solver = getSolver(args);
        if(args.size()>0) {
           solver.assertPB(args, c, compareTo);
        }else{
            int trueCount = 0;
            for(Lit l:args){
                if(l==Lit.True){
                    trueCount++;
                }
            }
            //statically check if the comparison holds
            if(!checkStaticPB(trueCount,c,compareTo)){
                contradiction();
            }
        }
    }
    public static void assertPB(List<Lit> args, List<Integer> weights, Comparison c, int compareTo) {
        Solver solver = getSolver(args);
        if(args.size()>0) {
            solver.assertPB(args, c, compareTo);
        }else{
            int trueCount = 0;
            for(int i = 0;i<args.size();i++){
                Lit l = args.get(i);
                int weight = 1;
                if(i<weights.size()){
                    weight = weights.get(i);
                }
                if(l==Lit.True){
                    trueCount+=weight;
                }
            }
            //statically check if the comparison holds
            if(!checkStaticPB(trueCount,c,compareTo)){
                contradiction();
            }
        }
    }

}
