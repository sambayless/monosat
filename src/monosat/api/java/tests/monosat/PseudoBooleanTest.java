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
import java.util.List;

import static monosat.Lit.False;
import static monosat.Lit.True;
import static org.junit.Assert.*;

/**
 * Tests psuedo-Boolean constraints
 */
public class PseudoBooleanTest {


  @Test
  public void test_GEQ() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    solver.assertPB(Arrays.asList(a,b,c), Comparison.GEQ,2);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.GEQ,3);
    assertTrue(solver.solve());
    solver.assertNand(a,b,c);
    assertFalse(solver.solve());
  }

  @Test
  public void test_GT() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    solver.assertPB(Arrays.asList(a,b,c), Comparison.GT,1);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.GT,2);
    assertTrue(solver.solve());
    solver.assertNand(a,b,c);
    assertFalse(solver.solve());
  }


  @Test
  public void test_GTEdgeCase() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    solver.assertPB(Arrays.asList(a,b,c), Comparison.GT,2);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.GT,3);
    assertFalse(solver.solve());
  }

  @Test
  public void test_LEQ() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LEQ,2);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LEQ,1);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LEQ,0);
    assertTrue(solver.solve());
    solver.assertOr(a,b,c);
    assertFalse(solver.solve());
  }


  @Test
  public void test_LT() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LT,3);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LT,2);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LT,1);
    assertTrue(solver.solve());
    solver.assertOr(a,b,c);
    assertFalse(solver.solve());
  }

  @Test
  public void test_LTEdgeCase() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LT,3);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LT,2);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LT,1);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.LT,0);
    assertFalse(solver.solve());
  }


  @Test
  public void test_EQ() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    solver.assertPB(Arrays.asList(a,b,c), Comparison.EQ,3);
    assertTrue(solver.solve());
    solver.assertPB(Arrays.asList(a,b,c), Comparison.EQ,2);
    assertFalse(solver.solve());
  }


  @Test
  public void test_NEQ() {
    Solver solver = new Solver();
    Lit a = new Lit(solver);
    Lit b = new Lit(solver);
    Lit c = new Lit(solver);
    try{
      solver.assertPB(Arrays.asList(a,b,c), Comparison.NEQ,3);
      assertTrue(solver.solve());
      solver.assertPB(Arrays.asList(a,b,c), Comparison.NEQ,2);
      assertTrue(solver.solve());
      solver.assertPB(Arrays.asList(a,b,c), Comparison.NEQ,1);
      assertTrue(solver.solve());
      solver.assertPB(Arrays.asList(a,b,c), Comparison.NEQ,0);
      assertFalse(solver.solve());
      fail("Expected UnsupportedOperationException");
    }catch(UnsupportedOperationException e){
      //NEQ not yet implemented for PB constraints
    }
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


}
