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

import org.junit.Test;

import java.util.ArrayList;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

/**
 * Tests the basic API functionality of the solver.
 * The goal of these tests is simply to make sure the JNI implementation works as expected,
 * but does not attempt to fully test the correctness of the underlying solver.
 */
public class LogicTest {

    @Test
    public void test_newLit() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        assertNotEquals(a, solver.True());

        Lit b = solver.newLit();
        assertNotEquals(a, b);
        assertNotEquals(a.l, b.l);
        assertNotEquals(a.toVar(), b.toVar());
        assertEquals(solver.solve(a), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(Logic.not(b)), true);
        assertEquals(solver.solve(a), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(a), b), true);
        assertEquals(solver.solve(Logic.not(a), a), false);
        assertEquals(solver.solve(a), true);
        assertEquals(solver.solve(b), true);
    }

    @Test
    public void test_solve() {
        Solver solver = new Solver();
        assertEquals(solver.solve(), true);
        assertEquals(solver.solve(), true);
        Logic.assertTrue(solver.True());
        assertEquals(solver.solve(), true);
        Logic.assertFalse(solver.True());
        assertEquals(solver.solve(), false);
        Solver solver2 = new Solver();
        assertEquals(solver2.solve(), true);
    }


    @Test
    public void test_testSolveAssumps() {
        Solver solver = new Solver();
        assertEquals(solver.solve(), true);
        assertEquals(solver.solve(), true);
        assertEquals(solver.solve(solver.True()), true);
        assertEquals(solver.solve(solver.False()), false);
        assertEquals(solver.solve(solver.True(), solver.False()), false);
        assertEquals(solver.solve(), true);
    }

    @Test
    public void test_ite() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit result = Logic.ite(c, a, b);
        assertEquals(solver.solve(c, a, Logic.not(b), result), true);
        assertEquals(solver.solve(c, a, Logic.not(b), Logic.not(result)), false);
        assertEquals(solver.solve(c, Logic.not(a), Logic.not(b), Logic.not(result)), true);

        assertEquals(solver.solve(Logic.not(c), a, Logic.not(b), result), false);
        assertEquals(solver.solve(Logic.not(c), a, Logic.not(b), Logic.not(result)), true);
        assertEquals(solver.solve(Logic.not(c), Logic.not(a), b, result), true);

    }

    @Test
    public void test_not() {
        Solver solver = new Solver();
        assertEquals(Logic.not(solver.True()), solver.False());
        Lit a = solver.newLit();
        assertNotEquals(Logic.not(a), a);
        assertEquals(Logic.not(a), a.negate());
    }

    @Test
    public void test_and() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit result = Logic.and(a, b);
        assertEquals(solver.solve(a, b, result), true);
        assertEquals(solver.solve(a, Logic.not(b), result), false);
        assertEquals(solver.solve(a, Logic.not(b), Logic.not(result)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)), true);
        assertEquals(solver.solve(a, b, result), true);
    }

    @Test
    public void test_ands() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit result = Logic.and(a, b, c);
        assertEquals(solver.solve(a, b, c, result), true);
        assertEquals(solver.solve(a, Logic.not(b), c, result), false);
        assertEquals(solver.solve(a, Logic.not(b), c, Logic.not(result)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)), true);
        assertEquals(solver.solve(a, b, c, result), true);
    }

    @Test
    public void test_or() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit result = Logic.or(a, b);
        assertEquals(solver.solve(a, b, result), true);
        assertEquals(solver.solve(a, Logic.not(b), result), true);
        assertEquals(solver.solve(a, Logic.not(b), Logic.not(result)), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), result), false);
        assertEquals(solver.solve(a, b, result), true);
    }

    @Test
    public void test_ors() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit result = Logic.or(a, b, c);
        assertEquals(solver.solve(a, b, c, result), true);
        assertEquals(solver.solve(a, Logic.not(b), c, result), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c), result), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c, result), true);
        assertEquals(solver.solve(a, b, c, result), true);
    }


    @Test
    public void test_nand() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit result = Logic.nand(a, b);
        assertEquals(solver.solve(a, b, result), false);
        assertEquals(solver.solve(a, Logic.not(b), result), true);
        assertEquals(solver.solve(a, Logic.not(b), Logic.not(result)), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)), false);
        assertEquals(solver.solve(a, b, result), false);
    }

    @Test
    public void test_nands() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit result = Logic.nand(a, b, c);
        assertEquals(solver.solve(a, b, c, result), false);
        assertEquals(solver.solve(a, Logic.not(b), c, result), true);
        assertEquals(solver.solve(a, Logic.not(b), c, Logic.not(result)), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)), false);
        assertEquals(solver.solve(a, b, c, result), false);
    }

    @Test
    public void test_nor() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit result = Logic.nor(a, b);
        assertEquals(solver.solve(a, b, result), false);
        assertEquals(solver.solve(a, Logic.not(b), result), false);
        assertEquals(solver.solve(a, Logic.not(b), Logic.not(result)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), result), true);
        assertEquals(solver.solve(a, b, result), false);
    }

    @Test
    public void test_nors() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit result = Logic.nor(a, b, c);
        assertEquals(solver.solve(a, b, c, result), false);
        assertEquals(solver.solve(a, Logic.not(b), c, result), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c), result), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c, result), false);
        assertEquals(solver.solve(a, b, c, result), false);
    }

    @Test
    public void test_xor() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit result = Logic.xor(a, b);
        assertEquals(solver.solve(a, b, result), false);
        assertEquals(solver.solve(a, Logic.not(b), result), true);
        assertEquals(solver.solve(a, Logic.not(b), Logic.not(result)), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)), true);
        assertEquals(solver.solve(a, b, result), false);
    }

    @Test
    public void test_xors() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit result = Logic.xor(a, b, c);
        assertEquals(solver.solve(a, b, c, result), true);
        assertEquals(solver.solve(a, Logic.not(b), c, result), false);
        assertEquals(solver.solve(a, Logic.not(b), c, Logic.not(result)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)), false);
        assertEquals(solver.solve(a, b, c, result), true);
    }

    @Test
    public void test_xnor() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit result = Logic.xnor(a, b);
        assertEquals(solver.solve(a, b, result), true);
        assertEquals(solver.solve(a, Logic.not(b), result), false);
        assertEquals(solver.solve(a, Logic.not(b), Logic.not(result)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)), false);
        assertEquals(solver.solve(a, b, result), true);
    }

    @Test
    public void test_xnors() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit result = Logic.xnor(a, b, c);
        assertEquals(solver.solve(a, b, c, result), false);
        assertEquals(solver.solve(a, Logic.not(b), c, result), true);
        assertEquals(solver.solve(a, Logic.not(b), c, Logic.not(result)), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)), true);
        assertEquals(solver.solve(a, b, c, result), false);
    }

    @Test
    public void test_assertTrue() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Logic.assertTrue(a);
        assertEquals(solver.solve(a), true);
        assertEquals(solver.solve(Logic.not(a)), false);
    }

    @Test
    public void test_assertFalse() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Logic.assertFalse(a);
        assertEquals(solver.solve(a), false);
        assertEquals(solver.solve(Logic.not(a)), true);
    }

    @Test
    public void test_assertAnd() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Logic.assertAnd(a, b);
        assertEquals(solver.solve(a, b), true);
        assertEquals(solver.solve(Logic.not(a)), false);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(b)), false);
    }

    @Test
    public void test_assertAnds() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Logic.assertAnd(a, b, c);
        assertEquals(solver.solve(a, b, c), true);
        assertEquals(solver.solve(Logic.not(a)), false);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(c)), false);
    }

    @Test
    public void test_assertOr() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Logic.assertOr(a, b);
        assertEquals(solver.solve(a, b), true);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(b)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b)), false);

    }

    @Test
    public void test_assertOrs() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Logic.assertOr(a, b, c);
        assertEquals(solver.solve(a, b, c), true);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(c)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c)), false);
    }

    @Test
    public void test_assertNand() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Logic.assertNand(a, b);
        assertEquals(solver.solve(a, b), false);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(b)), true);
    }

    @Test
    public void test_assertNands() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Logic.assertNand(a, b, c);
        assertEquals(solver.solve(a, b, c), false);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(c)), true);
    }

    @Test
    public void test_assertNor() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Logic.assertNor(a, b);
        assertEquals(solver.solve(a, b), false);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), false);
        assertEquals(solver.solve(Logic.not(b)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b)), true);

    }

    @Test
    public void test_assertNors() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Logic.assertNor(a, b, c);
        assertEquals(solver.solve(a, b, c), false);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), false);
        assertEquals(solver.solve(Logic.not(c)), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c)), true);
    }

    @Test
    public void test_assertXor() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Logic.assertXor(a, b);
        assertEquals(solver.solve(a, b), false);
        assertEquals(solver.solve(Logic.not(a), b), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b)), false);
        assertEquals(solver.solve(a, Logic.not(b)), true);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(b)), true);
    }

    @Test
    public void test_assertXors() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Logic.assertXor(a, b);
        assertEquals(solver.solve(a, b, c), false);
        assertEquals(solver.solve(Logic.not(a), b, c), true);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c), false);
        assertEquals(solver.solve(a, Logic.not(b), c), true);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(c), true);
        assertEquals(solver.solve(Logic.not(c)), true);
    }


    @Test
    public void test_assertXnor() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Logic.assertXnor(a, b);
        assertEquals(solver.solve(a, b), true);
        assertEquals(solver.solve(Logic.not(a), b), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b)), true);
        assertEquals(solver.solve(a, Logic.not(b)), false);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(b)), true);
    }

    @Test
    public void test_assertXnors() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Logic.assertXnor(a, b);
        assertEquals(solver.solve(a, b, c), true);
        assertEquals(solver.solve(Logic.not(a), b, c), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b), c), true);
        assertEquals(solver.solve(a, Logic.not(b), c), false);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(c), true);
        assertEquals(solver.solve(Logic.not(c)), true);
    }

    @Test
    public void test_assertEqual() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Logic.assertEqual(a, b);
        assertEquals(solver.solve(a, b), true);
        assertEquals(solver.solve(Logic.not(a), b), false);
        assertEquals(solver.solve(Logic.not(a), Logic.not(b)), true);
        assertEquals(solver.solve(a, Logic.not(b)), false);
        assertEquals(solver.solve(Logic.not(a)), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(Logic.not(b)), true);
    }

    @Test
    public void test_assertAtMostOne() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Logic.assertAtMostOne(a,b,c);

        assertEquals(solver.solve(a, b), false);
        assertEquals(solver.solve(b, c), false);
        assertEquals(solver.solve(a, c), false);
        assertEquals(solver.solve(a, b.negate()), true);
        assertEquals(solver.solve(b), true);
        assertEquals(solver.solve(c),true);
    }


    @Test
    public void test_unsatCore() {
        Solver solver = new Solver();
        Lit a = solver.newLit();
        Lit b = solver.newLit();
        Lit c = solver.newLit();
        Lit d = solver.newLit();
        Lit x = Logic.equal(a, b);
        Lit y = Logic.equal(b, c.negate());
        Lit z = Logic.equal(c, a);
        Lit w = Logic.equal(b, d.negate());
        Lit q = Logic.equal(d, a);
        assertEquals(solver.solve(), true);
        assertEquals(solver.solve(x, y,z,w,q), false);
        ArrayList<Lit> conflict = solver.getConflictClause();
        assertEquals(conflict.isEmpty(),false);
        ArrayList<Lit> test = new ArrayList<>();
        for(Lit l:conflict){
            assert(!l.sign());
            test.add(l.negate());
        }
        assertEquals(solver.solve(test), false);

        assertEquals(solver.solve(x, y,z,w,q), false);
        ArrayList<Lit> conflict2 = solver.getConflictClause(true);
        assertEquals(conflict2.isEmpty(),false);
        assertEquals(conflict2.size()<=conflict.size(), true);
        ArrayList<Lit> test2 = new ArrayList<>();
        for(Lit l:conflict2){
            assert(!l.sign());
            test2.add(l.negate());
        }
        assertEquals(solver.solve(test2), false);
        assertEquals(solver.solve(), true);

        ArrayList<Lit> conflict3 = solver.minimizeUnsatCore(x, y,z,w,q);
        assertEquals(conflict3.isEmpty(),false);
        assertEquals(conflict3.size()<=3, true);
        assertEquals(solver.solve(), true);
    }


}
