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

import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static monosat.Lit.False;
import static monosat.Lit.True;
import static org.junit.Assert.*;

/**
 * Tests the basic API functionality of the solver. The goal of these tests is simply to make sure
 * the JNI implementation works as expected, but does not attempt to fully test the correctness of
 * the underlying solver.
 */
public class LogicTest {

  @Test
  public void test_newLit() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    assertNotEquals(a, False);

    Lit b = new Lit(solver);
    assertNotEquals(a, b);
    assertNotEquals(a.l, b.l);
    assertNotEquals(a.toVar(), b.toVar());
    assertTrue(solver.solve(a));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(Logic.not(b)));
    assertTrue(solver.solve(a));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(a), b));
    assertFalse(solver.solve(Logic.not(a), a));
    assertTrue(solver.solve(a));
    assertTrue(solver.solve(b));
  }

  @Test
  public void test_solve() {
    Solver solver = new Solver();
    assertTrue(solver.solve());
    assertTrue(solver.solve());
    Logic.assertTrue(True);
    assertTrue(solver.solve());
    Logic.disallowContradictions();
    try {
      Logic.assertFalse(True);
      fail("Expected a contradiction exception");
    } catch (TrivialContradictionException e) {
      // ok
    }
    assertFalse(solver.solve());
    Solver solver2 = new Solver();
    assertTrue(solver2.solve());
  }

  @Test
  public void test_testSolveAssumps() {
    Solver solver = new Solver();
    assertTrue(solver.solve());
    assertTrue(solver.solve());
    assertTrue(solver.solve(True));
    assertFalse(solver.solve(False));
    assertFalse(solver.solve(True, False));
    assertTrue(solver.solve());
  }

  @Test
  public void test_ite() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.ite(c, a, b);
    assertTrue(solver.solve(c, a, Logic.not(b), result));
    assertFalse(solver.solve(c, a, Logic.not(b), Logic.not(result)));
    assertTrue(solver.solve(c, Logic.not(a), Logic.not(b), Logic.not(result)));

    assertFalse(solver.solve(Logic.not(c), a, Logic.not(b), result));
    assertTrue(solver.solve(Logic.not(c), a, Logic.not(b), Logic.not(result)));
    assertTrue(solver.solve(Logic.not(c), Logic.not(a), b, result));
  }

  @Test
  public void test_not() {
    Solver solver = new Solver();
    assertEquals(Logic.not(True), False);
    Lit a = new Lit(solver);
    assertNotEquals(Logic.not(a), a);
    assertEquals(Logic.not(a), a.not());
  }

  @Test
  public void test_and() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit result = Logic.and(a, b);
    assertTrue(solver.solve(a, b, result));
    assertFalse(solver.solve(a, Logic.not(b), result));
    assertTrue(solver.solve(a, Logic.not(b), Logic.not(result)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)));
    assertTrue(solver.solve(a, b, result));
  }

  @Test
  public void test_ands() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.and(a, b, c);
    assertTrue(solver.solve(a, b, c, result));
    assertFalse(solver.solve(a, Logic.not(b), c, result));
    assertTrue(solver.solve(a, Logic.not(b), c, Logic.not(result)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)));
    assertTrue(solver.solve(a, b, c, result));
  }

  @Test
  public void test_or() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit result = Logic.or(a, b);
    assertTrue(solver.solve(a, b, result));
    assertTrue(solver.solve(a, Logic.not(b), result));
    assertFalse(solver.solve(a, Logic.not(b), Logic.not(result)));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), result));
    assertTrue(solver.solve(a, b, result));
  }

  @Test
  public void test_ors() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.or(a, b, c);
    assertTrue(solver.solve(a, b, c, result));
    assertTrue(solver.solve(a, Logic.not(b), c, result));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c), result));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), c, result));
    assertTrue(solver.solve(a, b, c, result));
  }

  @Test
  public void test_nand() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit result = Logic.nand(a, b);
    assertFalse(solver.solve(a, b, result));
    assertTrue(solver.solve(a, Logic.not(b), result));
    assertFalse(solver.solve(a, Logic.not(b), Logic.not(result)));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)));
    assertFalse(solver.solve(a, b, result));
  }

  @Test
  public void test_nands() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.nand(a, b, c);
    assertFalse(solver.solve(a, b, c, result));
    assertTrue(solver.solve(a, Logic.not(b), c, result));
    assertFalse(solver.solve(a, Logic.not(b), c, Logic.not(result)));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)));
    assertFalse(solver.solve(a, b, c, result));
  }

  @Test
  public void test_nor() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit result = Logic.nor(a, b);
    assertFalse(solver.solve(a, b, result));
    assertFalse(solver.solve(a, Logic.not(b), result));
    assertTrue(solver.solve(a, Logic.not(b), Logic.not(result)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), result));
    assertFalse(solver.solve(a, b, result));
  }

  @Test
  public void test_nors() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.nor(a, b, c);
    assertFalse(solver.solve(a, b, c, result));
    assertFalse(solver.solve(a, Logic.not(b), c, result));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c), result));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), c, result));
    assertFalse(solver.solve(a, b, c, result));
  }

  @Test
  public void test_xor() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit result = Logic.xor(a, b);
    assertFalse(solver.solve(a, b, result));
    assertTrue(solver.solve(a, Logic.not(b), result));
    assertFalse(solver.solve(a, Logic.not(b), Logic.not(result)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)));
    assertFalse(solver.solve(a, b, result));
  }

  @Test
  public void test_xors() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.xor(a, b, c);
    assertTrue(solver.solve(a, b, c, result));
    assertFalse(solver.solve(a, Logic.not(b), c, result));
    assertTrue(solver.solve(a, Logic.not(b), c, Logic.not(result)));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)));
    assertTrue(solver.solve(a, b, c, result));
  }

  @Test
  public void test_xnor() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit result = Logic.xnor(a, b);
    assertTrue(solver.solve(a, b, result));
    assertFalse(solver.solve(a, Logic.not(b), result));
    assertTrue(solver.solve(a, Logic.not(b), Logic.not(result)));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), Logic.not(result)));
    assertTrue(solver.solve(a, b, result));
  }

  @Test
  public void test_xnors() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.xnor(a, b, c);
    assertFalse(solver.solve(a, b, c, result));
    assertTrue(solver.solve(a, Logic.not(b), c, result));
    assertFalse(solver.solve(a, Logic.not(b), c, Logic.not(result)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), c, Logic.not(result)));
    assertFalse(solver.solve(a, b, c, result));
  }

  @Test
  public void test_assertTrue() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Logic.assertTrue(a);
    assertTrue(solver.solve(a));
    assertFalse(solver.solve(Logic.not(a)));
  }

  @Test
  public void test_assertFalse() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Logic.assertFalse(a);
    assertFalse(solver.solve(a));
    assertTrue(solver.solve(Logic.not(a)));
  }

  @Test
  public void test_assertAnd() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Logic.assertAnd(a, b);
    assertTrue(solver.solve(a, b));
    assertFalse(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertFalse(solver.solve(Logic.not(b)));
  }

  @Test
  public void test_assertAnds() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Logic.assertAnd(a, b, c);
    assertTrue(solver.solve(a, b, c));
    assertFalse(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertFalse(solver.solve(Logic.not(c)));
  }

  @Test
  public void test_assertOr() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Logic.assertOr(a, b);
    assertTrue(solver.solve(a, b));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(b)));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b)));
  }

  @Test
  public void test_assertOrs() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Logic.assertOr(a, b, c);
    assertTrue(solver.solve(a, b, c));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(c)));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c)));
  }

  @Test
  public void test_assertNand() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Logic.assertNand(a, b);
    assertFalse(solver.solve(a, b));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(b)));
  }

  @Test
  public void test_assertNands() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Logic.assertNand(a, b, c);
    assertFalse(solver.solve(a, b, c));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(c)));
  }

  @Test
  public void test_assertNor() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Logic.assertNor(a, b);
    assertFalse(solver.solve(a, b));
    assertTrue(solver.solve(Logic.not(a)));
    assertFalse(solver.solve(b));
    assertTrue(solver.solve(Logic.not(b)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b)));
  }

  @Test
  public void test_assertNors() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Logic.assertNor(a, b, c);
    assertFalse(solver.solve(a, b, c));
    assertTrue(solver.solve(Logic.not(a)));
    assertFalse(solver.solve(b));
    assertTrue(solver.solve(Logic.not(c)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), Logic.not(c)));
  }

  @Test
  public void test_assertXor() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Logic.assertXor(a, b);
    assertFalse(solver.solve(a, b));
    assertTrue(solver.solve(Logic.not(a), b));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b)));
    assertTrue(solver.solve(a, Logic.not(b)));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(b)));
  }

  @Test
  public void test_assertXors() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Logic.assertXor(a, b);
    assertFalse(solver.solve(a, b, c));
    assertTrue(solver.solve(Logic.not(a), b, c));
    assertFalse(solver.solve(Logic.not(a), Logic.not(b), c));
    assertTrue(solver.solve(a, Logic.not(b), c));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(c));
    assertTrue(solver.solve(Logic.not(c)));
  }

  @Test
  public void test_assertXnor() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Logic.assertXnor(a, b);
    assertTrue(solver.solve(a, b));
    assertFalse(solver.solve(Logic.not(a), b));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b)));
    assertFalse(solver.solve(a, Logic.not(b)));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(b)));
  }

  @Test
  public void test_assertXnors() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Logic.assertXnor(a, b);
    assertTrue(solver.solve(a, b, c));
    assertFalse(solver.solve(Logic.not(a), b, c));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), c));
    assertFalse(solver.solve(a, Logic.not(b), c));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(c));
    assertTrue(solver.solve(Logic.not(c)));
  }

  @Test
  public void test_assertEqual() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Logic.assertEqual(a, b);
    assertTrue(solver.solve(a, b));
    assertFalse(solver.solve(Logic.not(a), b));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b)));
    assertFalse(solver.solve(a, Logic.not(b)));
    assertTrue(solver.solve(Logic.not(a)));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(Logic.not(b)));
  }

  @Test
  public void test_assertAtMostOne() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Logic.assertAtMostOne(a, b, c);

    assertFalse(solver.solve(a, b));
    assertFalse(solver.solve(b, c));
    assertFalse(solver.solve(a, c));
    assertTrue(solver.solve(a, b.not()));
    assertTrue(solver.solve(b));
    assertTrue(solver.solve(c));

    Logic.assertAtMostOne();//this should have no effect
    assertTrue(solver.solve(c));

    solver.assertAtMostOne();
    assertTrue(solver.solve(c));
  }

  @Test
  public void test_unsatCore() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit d = new Lit(solver);
    Lit x = Logic.equal(a, b);
    Lit y = Logic.equal(b, c.not());
    Lit z = Logic.equal(c, a);
    Lit w = Logic.equal(b, d.not());
    Lit q = Logic.equal(d, a);
    assertTrue(solver.solve());
    assertFalse(solver.solve(x, y, z, w, q));
    List<Lit> conflict = solver.getConflictClause();
    assertFalse(conflict.isEmpty());
    ArrayList<Lit> test = new ArrayList<>();
    for (Lit l : conflict) {
      assert (!l.sign());
      test.add(l.not());
    }
    assertFalse(solver.solve(test));

    assertFalse(solver.solve(x, y, z, w, q));
    List<Lit> conflict2 = solver.getConflictClause(true);
    assertFalse(conflict2.isEmpty());
    assertTrue(conflict2.size() <= conflict.size());
    ArrayList<Lit> test2 = new ArrayList<>();
    for (Lit l : conflict2) {
      assert (!l.sign());
      test2.add(l.not());
    }
    assertFalse(solver.solve(test2));
    assertTrue(solver.solve());

    List<Lit> conflict3 = solver.minimizeUnsatCore(x, y, z, w, q);
    assertFalse(conflict3.isEmpty());
    assertTrue(conflict3.size() <= 3);
    assertTrue(solver.solve());
  }

  @Test
  public void test_Constant() {
    Solver solver = new Solver();
    assertTrue(True.isConst());
    assertTrue(True.isConstTrue());
    assertFalse(True.isConstFalse());

    assertTrue(False.isConst());
    assertFalse(False.isConstTrue());
    assertTrue(False.isConstFalse());

    Lit a = new Lit(solver);

    assertFalse(a.isConst());
    assertFalse(a.isConstTrue());
    assertFalse(a.isConstFalse());

    Lit b = new Lit(solver);
    solver.assertTrue(a);
    assertTrue(a.isConst());
    assertTrue(a.isConstTrue());
    assertFalse(a.isConstFalse());

    assertFalse(b.isConst());
    assertFalse(b.isConstTrue());
    assertFalse(b.isConstFalse());
  }

  @Test
  public void test_AMO() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    Lit result = Logic.xnor(a, b, c);
    assertTrue(solver.solve());
    // An empty AMO call should have no effect
    Logic.assertAtMostOne();
    assertTrue(solver.solve());
    Logic.assertAtMostOne(a);
    assertTrue(solver.solve());
    assertTrue(solver.solve(a));
    assertTrue(solver.solve(a.not()));

    Logic.assertAtMostOne(a, b);
    assertTrue(solver.solve());
    assertTrue(solver.solve(a));
    assertTrue(solver.solve(a.not(), b.not()));
  }

    @Test
    public void testLoading() throws IOException {
        File f = File.createTempFile("test", ".gnf");
        String filename = f.getAbsolutePath().toString();
        f.delete();

        {
          Solver solver = new monosat.Solver("", filename);

          Lit a = new Lit(solver,"a");
          Lit b = new Lit(solver,"b");
          Lit c = new Lit(solver,"c");
          Logic.assertXnor(a, b);
          assertTrue(solver.solve(a, b, c));
          assertFalse(solver.solve(Logic.not(a), b, c));
          assertTrue(solver.solve(Logic.not(a), Logic.not(b), c));
          assertFalse(solver.solve(a, Logic.not(b), c));
          assertTrue(solver.solve(Logic.not(a)));
          assertTrue(solver.solve(c));
          assertTrue(solver.solve(Logic.not(c)));
        }
        monosat.Solver solver = new monosat.Solver();
        solver.loadConstraints(filename);

        Lit a = solver.getLiteral("a");
        Lit b = solver.getLiteral("b");
        Lit c = solver.getLiteral("c");

        assertTrue(solver.solve());
        assertTrue(solver.solve(a, b, c));
        assertFalse(solver.solve(Logic.not(a), b, c));
        assertTrue(solver.solve(Logic.not(a), Logic.not(b), c));
        assertFalse(solver.solve(a, Logic.not(b), c));
        assertTrue(solver.solve(Logic.not(a)));
        assertTrue(solver.solve(c));
        assertTrue(solver.solve(Logic.not(c)));
    }


  @Test
  public void testLoadingAMO() throws IOException {
    File f = File.createTempFile("test", ".gnf");
    String filename = f.getAbsolutePath().toString();
    f.delete();

    {
      Solver solver = new monosat.Solver("", filename);

      Lit a = new Lit(solver,"a");
      Lit b = new Lit(solver,"b");
      Lit c = new Lit(solver,"c");
      Logic.assertAtMostOne(a, b,c);
      assertFalse(solver.solve(a, b, c));
      assertTrue(solver.solve(a, Logic.not(b), c.not()));
      assertTrue(solver.solve(Logic.not(a), b, Logic.not(c)));
      assertTrue(solver.solve(Logic.not(a), Logic.not(b), c));
    }
    monosat.Solver solver = new monosat.Solver();
    solver.loadConstraints(filename);

    Lit a = solver.getLiteral("a");
    Lit b = solver.getLiteral("b");
    Lit c = solver.getLiteral("c");

    assertTrue(solver.solve());
    assertFalse(solver.solve(a, b, c));
    assertTrue(solver.solve(a, Logic.not(b), c.not()));
    assertTrue(solver.solve(Logic.not(a), b, Logic.not(c)));
    assertTrue(solver.solve(Logic.not(a), Logic.not(b), c));
  }
}
