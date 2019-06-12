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
import java.util.Arrays;

import static org.junit.Assert.*;

public class SolverTest {
  @Test
  public void testSolve() {
    monosat.Solver s = new monosat.Solver();
    assertTrue(s.solve());
    assertTrue(s.solve());
    assertTrue(s.solve(monosat.Lit.True));
    assertFalse(s.solve(monosat.Lit.False));
    assertFalse(s.solve(monosat.Lit.True, monosat.Lit.False));
    assertTrue(s.solve());
  }

  @Test
  public void testConstraints() {
    monosat.Solver s = new monosat.Solver();
    Constraints.nqueens(s,4);
    assertTrue(s.solve());

  }

  @Test
  public void testTryWithResources() {
    try(monosat.Solver s = new monosat.Solver()) {
      assertTrue(s.solve());
    }
  }

  @Test
  public void testArguments() {
    monosat.Solver s = new monosat.Solver("-no-reach-underapprox-cnf");
    assertTrue(s.solve());
    assertTrue(s.solve());
    assertTrue(s.solve(monosat.Lit.True));
    assertFalse(s.solve(monosat.Lit.False));
    assertFalse(s.solve(monosat.Lit.True, monosat.Lit.False));
    assertTrue(s.solve());
  }

  @Test
  public void testMultipleSolvers() {
    monosat.Solver s = new monosat.Solver();
    assertTrue(s.solve());

    monosat.Solver s2 = new monosat.Solver();
    assertTrue(s2.solve());
    assertTrue(s.solve());
    s.close();

    assertTrue(s2.solve());
  }

  @Test
  public void testLoad() throws IOException {
    File f = File.createTempFile("test", ".gnf");
    String filename = f.getAbsolutePath().toString();
    f.delete();
    monosat.Solver s = new monosat.Solver("", filename);
    Lit a = new Lit(s);
    s.addClause(a.not());
    assertTrue(s.solve());
    assertFalse(s.solve(a));
    assertTrue(s.solve(a.not()));
    s.flushConstraintFile();//these are both optional
    s.closeConstraintFile();//these are both optional
    s.close();

    monosat.Solver s2 = new monosat.Solver();
    assertTrue(s2.solve());
    s2.loadConstraints(filename);
    assertTrue(s2.solve());
  }

  @Test
  public void testOrUnsat() {
    monosat.Solver s = new monosat.Solver();
    assertTrue(s.solve());
    s.assertAnd();
    assertTrue(s.solve());
    s.assertOr();
    assertFalse(s.solve());
    monosat.Solver s2 = new monosat.Solver();
    assertTrue(s2.solve());
    assertFalse(s.solve());
  }

  @Test
  public void testAddClause() {
    monosat.Solver s = new monosat.Solver();
    assertTrue(s.solve());
    s.addClause();
    assertFalse(s.solve());
    s = new monosat.Solver();
    monosat.Lit a = new monosat.Lit(s);
    monosat.Lit b = new monosat.Lit(s);
    monosat.Lit c = new monosat.Lit(s);
    monosat.Lit d = new monosat.Lit(s);
    s.addClause(a);
    s.addClause(a, b.not());
    s.addClause(a.not(), b, c);
    s.addClause(a.not(), b, c, d.not());
    assertTrue(s.solve());
    ArrayList<monosat.Lit> clause =
        new ArrayList<monosat.Lit>(Arrays.asList(a, b.not(), c.not(), d));
    s.addClause(clause);
    assertTrue(s.solve());
  }

  @Test
  public void test_AMO() {
    monosat.Solver solver = new monosat.Solver();
    monosat.Lit a = new monosat.Lit(solver);
    monosat.Lit b = new monosat.Lit(solver);
    monosat.Lit c = new monosat.Lit(solver);
    monosat.Lit result = monosat.Logic.xnor(a, b, c);
    assertTrue(solver.solve());
    // An empty AMO call should have no effect
    solver.assertAtMostOne();
    assertTrue(solver.solve());
    solver.assertAtMostOne(a);
    assertTrue(solver.solve());
    assertTrue(solver.solve(a));
    assertTrue(solver.solve(a.not()));

    solver.assertAtMostOne(a, b);
    assertTrue(solver.solve());
    assertTrue(solver.solve(a));
    assertTrue(solver.solve(a.not(), b.not()));
  }

  @Test
  public void getLits() {
    monosat.Solver s = new monosat.Solver();
    assertEquals(s.getLits().size(), 1);
    monosat.Lit a = new monosat.Lit(s);
    monosat.Lit b = new monosat.Lit(s);
    monosat.Lit c = new monosat.Lit(s);
    monosat.Lit d = new monosat.Lit(s);
    assertEquals(s.getLits().size(), 5);
    s.addClause(a);
    s.addClause(a, b.not());
    s.addClause(a.not(), b, c);
    s.addClause(a.not(), b, c, d.not());
    assertTrue(s.solve());
    ArrayList<monosat.Lit> clause =
        new ArrayList<monosat.Lit>(Arrays.asList(a, b.not(), c.not(), d));
    s.addClause(clause);
    assertTrue(s.solve());

    assertEquals(s.getLits().size(), 5); // plus 1, for the constant true literal
    for (monosat.Lit l : s.getLits()) {
      boolean val = l.value();
    }
  }

  @Test
  public void testOk() {
    monosat.Solver s = new monosat.Solver();
    assertTrue(s.ok());
    assertTrue(s.solve());
    assertTrue(s.ok());
    assertTrue(s.solve(monosat.Lit.True));
    assertTrue(s.ok());
    assertFalse(s.solve(monosat.Lit.False));
    assertTrue(s.ok());
    s.addClause(monosat.Lit.False);
    assertFalse(s.ok());
    assertFalse(s.solve());
    assertFalse(s.ok());
    s.close();
    try{
      s.ok();
      fail("Accessing native solver methods after closing the solver should throw an NPE");
    }catch(NullPointerException e){
      //ok
    }
    //but creating a new solver should be ok
    monosat.Solver s2 = new monosat.Solver();
    assertTrue(s2.ok());
    assertTrue(s2.solve());
    assertTrue(s2.ok());
 }

  @Test
  public void version() {
    String version = monosat.Solver.getVersion();
    System.out.println(version);
    assertFalse(version.isEmpty());
  }
}
