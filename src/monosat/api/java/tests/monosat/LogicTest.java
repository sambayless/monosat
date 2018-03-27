package monosat;
import monosat.Lit;
import monosat.Solver;
import org.junit.Test;

import monosat.Logic;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Tests the basic API functionality of the solver.
 * The goal of these tests is simply to make sure the JNI implementation works as expected,
 * but does not attempt to fully test the correctness of the underlying solver.
 */
public class LogicTest {

    @Test
    public void test_newLit() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        assertNotEquals(a,Logic.True());

        Lit b = Logic.newLit();
        assertNotEquals(a,b);
        assertNotEquals(a.l,b.l);
        assertNotEquals(a.toVar(),b.toVar());
        assertEquals(Logic.solve(a),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(Logic.not(b)),true);
        assertEquals(Logic.solve(a),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(a),b),true);
        assertEquals(Logic.solve(Logic.not(a),a),false);
        assertEquals(Logic.solve(a),true);
        assertEquals(Logic.solve(b),true);
    }

    @Test
    public void test_solve() {
        Logic.newSolver();
        assertEquals(Logic.solve(), true);
        assertEquals(Logic.solve(), true);
        Logic.assertTrue(Logic.True());
        assertEquals(Logic.solve(), true);
        Logic.assertFalse(Logic.True());
        assertEquals(Logic.solve(), false);
        Logic.newSolver();
        assertEquals(Logic.solve(), true);
    }


    @Test
    public void test_testSolveAssumps() {
        Logic.newSolver();
        assertEquals(Logic.solve(), true);
        assertEquals(Logic.solve(), true);
        assertEquals(Logic.solve(Logic.True()), true);
        assertEquals(Logic.solve(Logic.False()), false);
        assertEquals(Logic.solve(Logic.True(), Logic.False()), false);
        assertEquals(Logic.solve(), true);
    }

    @Test
    public void test_ite() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Lit result = Logic.ite(c,a,b);
        assertEquals(Logic.solve(c, a,Logic.not(b),result),true);
        assertEquals(Logic.solve(c, a,Logic.not(b),Logic.not(result)),false);
        assertEquals(Logic.solve(c, Logic.not(a),Logic.not(b),Logic.not(result)),true);

        assertEquals(Logic.solve(Logic.not(c), a,Logic.not(b),result),false);
        assertEquals(Logic.solve(Logic.not(c), a,Logic.not(b),Logic.not(result)),true);
        assertEquals(Logic.solve(Logic.not(c), Logic.not(a),b,result),true);

    }

    @Test
    public void test_not() {
        assertEquals(Logic.not(Logic.True()),Logic.False());
        Lit a = Logic.newLit();
        assertNotEquals(Logic.not(a),a);
        assertEquals(Logic.not(a),a.negate());
    }

    @Test
    public void test_and() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit result = Logic.and(a,b);
        assertEquals(Logic.solve(a,b,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),result),false);
        assertEquals(Logic.solve(a,Logic.not(b),Logic.not(result)),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(result)),true);
        assertEquals(Logic.solve(a,b,result),true);
    }
    @Test
    public void test_ands() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Lit result = Logic.and(a,b,c);
        assertEquals(Logic.solve(a,b,c,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),c,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),c,Logic.not(result)),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c,Logic.not(result)),true);
        assertEquals(Logic.solve(a,b,c,result),true);
    }

    @Test
    public void test_or() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit result = Logic.or(a,b);
        assertEquals(Logic.solve(a,b,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),result),true);
        assertEquals(Logic.solve(a,Logic.not(b),Logic.not(result)),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),result),false);
        assertEquals(Logic.solve(a,b,result),true);
    }
    @Test
    public void test_ors() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Lit result = Logic.or(a,b,c);
        assertEquals(Logic.solve(a,b,c,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),c,result),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(c),result),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c,result),true);
        assertEquals(Logic.solve(a,b,c,result),true);
    }


    @Test
    public void test_nand() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit result = Logic.nand(a,b);
        assertEquals(Logic.solve(a,b,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),result),true);
        assertEquals(Logic.solve(a,Logic.not(b),Logic.not(result)),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(result)),false);
        assertEquals(Logic.solve(a,b,result),false);
    }
    @Test
    public void test_nands() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Lit result = Logic.nand(a,b,c);
        assertEquals(Logic.solve(a,b,c,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),c,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),c,Logic.not(result)),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c,Logic.not(result)),false);
        assertEquals(Logic.solve(a,b,c,result),false);
    }
    @Test
    public void test_nor() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit result = Logic.nor(a,b);
        assertEquals(Logic.solve(a,b,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),result),false);
        assertEquals(Logic.solve(a,Logic.not(b),Logic.not(result)),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),result),true);
        assertEquals(Logic.solve(a,b,result),false);
    }

    @Test
    public void test_nors() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Lit result = Logic.nor(a,b,c);
        assertEquals(Logic.solve(a,b,c,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),c,result),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(c),result),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c,result),false);
        assertEquals(Logic.solve(a,b,c,result),false);
    }
    @Test
    public void test_xor() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit result = Logic.xor(a,b);
        assertEquals(Logic.solve(a,b,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),result),true);
        assertEquals(Logic.solve(a,Logic.not(b),Logic.not(result)),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(result)),true);
        assertEquals(Logic.solve(a,b,result),false);
    }
    @Test
    public void test_xors() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Lit result = Logic.xor(a,b,c);
        assertEquals(Logic.solve(a,b,c,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),c,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),c,Logic.not(result)),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c,Logic.not(result)),false);
        assertEquals(Logic.solve(a,b,c,result),true);
    }
    @Test
    public void test_xnor() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit result = Logic.xnor(a,b);
        assertEquals(Logic.solve(a,b,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),result),false);
        assertEquals(Logic.solve(a,Logic.not(b),Logic.not(result)),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(result)),false);
        assertEquals(Logic.solve(a,b,result),true);
    }
    @Test
    public void test_xnors() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Lit result = Logic.xnor(a,b,c);
        assertEquals(Logic.solve(a,b,c,result),false);
        assertEquals(Logic.solve(a,Logic.not(b),c,result),true);
        assertEquals(Logic.solve(a,Logic.not(b),c,Logic.not(result)),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c,Logic.not(result)),true);
        assertEquals(Logic.solve(a,b,c,result),false);
    }
    @Test
    public void test_assertTrue() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Logic.assertTrue(a);
        assertEquals(Logic.solve(a),true);
        assertEquals(Logic.solve(Logic.not(a)),false);
    }

    @Test
    public void test_assertFalse() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Logic.assertFalse(a);
        assertEquals(Logic.solve(a),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
    }

    @Test
    public void test_assertAnd() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Logic.assertAnd(a,b);
        assertEquals(Logic.solve(a,b),true);
        assertEquals(Logic.solve(Logic.not(a)),false);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(b)),false);
    }
    @Test
    public void test_assertAnds() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Logic.assertAnd(a,b,c);
        assertEquals(Logic.solve(a,b,c),true);
        assertEquals(Logic.solve(Logic.not(a)),false);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(c)),false);
    }

    @Test
    public void test_assertOr() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Logic.assertOr(a,b);
        assertEquals(Logic.solve(a,b),true);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(b)),true);
        assertEquals(Logic.solve(Logic.not(a), Logic.not(b)),false);

    }
    @Test
    public void test_assertOrs() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Logic.assertOr(a,b,c);
        assertEquals(Logic.solve(a,b,c),true);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(c)),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(c)),false);
    }

    @Test
    public void test_assertNand() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Logic.assertNand(a,b);
        assertEquals(Logic.solve(a,b),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(b)),true);
    }
    @Test
    public void test_assertNands() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Logic.assertNand(a,b,c);
        assertEquals(Logic.solve(a,b,c),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(c)),true);
    }

    @Test
    public void test_assertNor() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Logic.assertNor(a,b);
        assertEquals(Logic.solve(a,b),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),false);
        assertEquals(Logic.solve(Logic.not(b)),true);
        assertEquals(Logic.solve(Logic.not(a), Logic.not(b)),true);

    }
    @Test
    public void test_assertNors() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Logic.assertNor(a,b,c);
        assertEquals(Logic.solve(a,b,c),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),false);
        assertEquals(Logic.solve(Logic.not(c)),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),Logic.not(c)),true);
    }

    @Test
    public void test_assertXor() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Logic.assertXor(a,b);
        assertEquals(Logic.solve(a,b),false);
        assertEquals(Logic.solve(Logic.not(a),b),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b)),false);
        assertEquals(Logic.solve(a,Logic.not(b)),true);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(b)),true);
    }
    @Test
    public void test_assertXors() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Logic.assertXor(a,b);
        assertEquals(Logic.solve(a,b,c),false);
        assertEquals(Logic.solve(Logic.not(a),b,c),true);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c),false);
        assertEquals(Logic.solve(a,Logic.not(b),c),true);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(c),true);
        assertEquals(Logic.solve(Logic.not(c)),true);
    }


    @Test
    public void test_assertXnor() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Logic.assertXnor(a,b);
        assertEquals(Logic.solve(a,b),true);
        assertEquals(Logic.solve(Logic.not(a),b),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b)),true);
        assertEquals(Logic.solve(a,Logic.not(b)),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(b)),true);
    }
    @Test
    public void test_assertXnors() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Lit c = Logic.newLit();
        Logic.assertXnor(a,b);
        assertEquals(Logic.solve(a,b,c),true);
        assertEquals(Logic.solve(Logic.not(a),b,c),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b),c),true);
        assertEquals(Logic.solve(a,Logic.not(b),c),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(c),true);
        assertEquals(Logic.solve(Logic.not(c)),true);
    }

    @Test
    public void test_assertEqual() {
        Logic.newSolver();
        Lit a = Logic.newLit();
        Lit b = Logic.newLit();
        Logic.assertEqual(a,b);
        assertEquals(Logic.solve(a,b),true);
        assertEquals(Logic.solve(Logic.not(a),b),false);
        assertEquals(Logic.solve(Logic.not(a),Logic.not(b)),true);
        assertEquals(Logic.solve(a,Logic.not(b)),false);
        assertEquals(Logic.solve(Logic.not(a)),true);
        assertEquals(Logic.solve(b),true);
        assertEquals(Logic.solve(Logic.not(b)),true);
    }



}
