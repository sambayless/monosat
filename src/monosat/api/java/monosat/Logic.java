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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

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
    //Literal level constructs
    public static Lit ite(Lit condition, Lit then, Lit els) {
        Solver solver = condition.getSolver();
        return solver.ite(condition, then, els);
    }

    public static Lit and(Lit... args) {
        if(args.length>=0){
            return args[0].getSolver().and(args);
        }else{
            throw new IllegalArgumentException("and requires at least one argument");
        }
    }

    public static Lit or(Lit... args) {
        if(args.length>=0){
            return args[0].getSolver().or(args);
        }else{
            throw new IllegalArgumentException("or requires at least one argument");
        }
    }

    public static Lit not(Lit a) {
        //don't need the solver for this call
        //so avoid the thread local access of using getSolver()
        return a.negate();
    }

    public static Lit nand(Lit... args) {
        if(args.length>=0){
            return args[0].getSolver().nand(args);
        }else{
            throw new IllegalArgumentException("nand requires at least one argument");
        }
    }

    public static Lit nor(Lit... args) {
        if(args.length>=0){
            return args[0].getSolver().nor(args);
        }else{
            throw new IllegalArgumentException("nor requires at least one argument");
        }
    }

    public static Lit xor(Lit... args) {
        if(args.length>=0){
            return args[0].getSolver().xor(args);
        }else{
            throw new IllegalArgumentException("xor requires at least one argument");
        }
    }

    public static Lit xnor(Lit... args) {
        if(args.length>0) {
            return args[0].getSolver().xnor(args);
        }else{
            throw new IllegalArgumentException("xnor requires at least one argument");
        }
    }

    public static Lit equal(Lit... args) {
        if(args.length>0) {
            return args[0].getSolver().xnor(args);
        }else{
            throw new IllegalArgumentException("equal requires at least one argument");
        }
    }
    //Assertion forms.
    public static void assertTrue(Lit a) {
        a.getSolver().assertTrue(a);
    }

    public static void assertFalse(Lit a) {
        a.getSolver().assertFalse(a);
    }

    public static void assertAnd(Lit... args) {
        if(args.length>0) {
            args[0].getSolver().assertAnd(args);
        }else{
            throw new IllegalArgumentException("assertAnd requires at least one argument");
        }
    }

    public static void assertOr(Lit... args) {
        if(args.length>0) {
            args[0].getSolver().assertOr(args);
        }else{
            throw new IllegalArgumentException("assertOr requires at least one argument");
        }
    }

    public static void assertNand(Lit... args) {
        if(args.length>0) {
            args[0].getSolver().assertNand(args);
        }else{
            throw new IllegalArgumentException("assertNand requires at least one argument");
        }
    }

    public static void assertNor(Lit... args) {
        if(args.length>0) {
            args[0].getSolver().assertNor(args);
        }else{
            throw new IllegalArgumentException("assertNor requires at least one argument");
        }
    }

    public static void assertXor(Lit... args) {
        if(args.length>0) {
            args[0].getSolver().assertXor(args);
        }else{
            throw new IllegalArgumentException("assertXor requires at least one argument");
        }
    }

    public static void assertXnor(Lit... args) {
        if(args.length>0) {
            args[0].getSolver().assertXnor(args);
        }else{
            throw new IllegalArgumentException("assertXnor requires at least one argument");
        }
    }

    public static void assertEqual(Lit a, Lit b) {
        a.getSolver().assertEqual(a, b);
    }

    public static void assertImplies(Lit a, Lit b) {a.getSolver().assertImplies(a,b);}

    //Multi-Literal constructs
    public static Lit and(Collection<Lit> args) {
        if(args.size()>0) {
            return args.iterator().next().getSolver().and(args);
        }else{
            throw new IllegalArgumentException("and requires at least one argument");
        }
    }

    public static Lit or(Collection<Lit> args) {
        if(args.size()>0) {
            return args.iterator().next().getSolver().or(args);
        }else{
            throw new IllegalArgumentException("or requires at least one argument");
        }
    }

    public static Lit nand(Collection<Lit> args) {
        if(args.size()>0) {
            return args.iterator().next().getSolver().nand(args);
        }else{
            throw new IllegalArgumentException("nand requires at least one argument");
        }
    }

    public static Lit nor(Collection<Lit> args) {
        if(args.size()>0) {
            return args.iterator().next().getSolver().nor(args);
        }else{
            throw new IllegalArgumentException("nor requires at least one argument");
        }
    }

    public static Lit xor(Collection<Lit> args) {
        if(args.size()>0) {
            return args.iterator().next().getSolver().xor(args);
        }else{
            throw new IllegalArgumentException("xor requires at least one argument");
        }
    }

    public static Lit xnor(Collection<Lit> args) {
        if(args.size()>0) {
            return args.iterator().next(). getSolver().xnor(args);
        }else{
            throw new IllegalArgumentException("xnor requires at least one argument");
        }
    }

    public static Lit implies(Lit a, Lit b){
        return a.getSolver().implies(a, b);
    }

    //assertion forms
    public static void assertAnd(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertAnd(args);
        }else{
            throw new IllegalArgumentException("assertAnd requires at least one argument");
        }
    }

    public static void assertOr(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertOr(args);
        }else{
            throw new IllegalArgumentException("assertOr requires at least one argument");
        }
    }

    public static void assertNand(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertNand(args);
        }else{
            throw new IllegalArgumentException("assertNand requires at least one argument");
        }
    }

    public static void assertNor(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertNor(args);
        }else{
            throw new IllegalArgumentException("assertNor requires at least one argument");
        }
    }

    public static void assertXor(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertXor(args);
        }else{
            throw new IllegalArgumentException("assertXor requires at least one argument");
        }
    }

    public static void assertXnor(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertXnor(args);
        }else{
            throw new IllegalArgumentException("assertXnor requires at least one argument");
        }
    }


    //Bitvector constructs
    public static BitVector ite(Lit condition, BitVector then, BitVector els) {
        return condition.getSolver().ite(condition, then, els);
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
        if(args.size()>0) {
            args.iterator().next().getSolver().maximizeLits(args);
        }else{
            throw new IllegalArgumentException("maximizeLits requires at least one argument");
        }
    }

    public static void minimizeLits(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().minimizeLits(args);
        }else{
            throw new IllegalArgumentException("minimizeLits requires at least one argument");
        }
    }

    public static void maximizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        if(literals.size()>0) {
            literals.iterator().next().getSolver().maximizeWeightedLits(literals, weights);
        }else{
            //do nothing
            //throw new IllegalArgumentException("maximizeWeightedLits requires at least one argument");
        }
    }

    public static void minimizeWeightedLits(Collection<Lit> literals, Collection<Integer> weights) {
        if(literals.size()>0) {
            literals.iterator().next().getSolver().minimizeWeightedLits(literals, weights);
        }else{
            //do nothing
            //throw new IllegalArgumentException("minimizeWeightedLits requires at least one argument");
        }
    }

    public static void assertAtMostOne(Collection<Lit> args) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertAtMostOne(args);
        }else{
            //do nothing
            //throw new IllegalArgumentException("assertAtMostOne requires at least one argument");
        }
    }

    public static void assertAtMostOne(Lit... args) {
        if(args.length>0) {
            args[0].getSolver().assertAtMostOne(args);
        }else{
            //do nothing
            //throw new IllegalArgumentException("assertAtMostOne requires at least one argument");
        }
    }
    public static void assertPB(Collection<Lit> args, Comparison c, int compareTo) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertPB(args, c, compareTo);
        }else{
            //do nothing
            //throw new IllegalArgumentException("assertPB requires at least one argument");
        }
    }
    public static void assertPB(Collection<Lit> args, Collection<Integer> weights, Comparison c, int compareTo) {
        if(args.size()>0) {
            args.iterator().next().getSolver().assertPB(args, weights, c, compareTo);
        }else{
            //do nothing
            //throw new IllegalArgumentException("assertPB requires at least one argument");
        }
    }

}
