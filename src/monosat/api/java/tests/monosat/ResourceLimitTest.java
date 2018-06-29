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

import org.junit.Assume;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Optional;

import static org.junit.Assert.*;

public class ResourceLimitTest {
  @Test
  public void testConflictLimit() {

    Solver s = new Solver();

    Constraints.unsatQueens(s,12);
    s.setConflictLimit(-1);
    assertFalse(s.solve());

    long n_conflicts =  s.nConflicts();
    assertTrue(n_conflicts>10);
    Solver s2 = new Solver();
    Constraints.unsatQueens(s2,12);
    s2.setConflictLimit(10);
    s2.setPropagationLimit(-1);
    assertEquals(s2.solveLimited(), Optional.empty());
    //limit still applies for a second solve
    assertEquals(s2.solveLimited(), Optional.empty());

    s2.setConflictLimit(-1);//disable limit and solve again
    assertFalse(s2.solveLimited().get());


  }

  @Test
  public void testSolveIgnoresLimits() {
    Solver s = new Solver();
    Constraints.unsatQueens(s,12);
    s.setConflictLimit(10);
    s.setPropagationLimit(10);
    //solver.solve() ignores resource limitations
    assertFalse(s.solve());
  }

  @Test
  public void testPropagationLimit() {

    Solver s = new Solver();

    Constraints.unsatQueens(s,12);
    s.setPropagationLimit(-1);
    assertFalse(s.solve());

    long n_props =  s.nPropagations();
    assertTrue(n_props>10);
    Solver s2 = new Solver();
    Constraints.unsatQueens(s2,12);
    s2.setPropagationLimit(10);
    s2.setConflictLimit(-1);
    assertEquals(s2.solveLimited(), Optional.empty());
    //limit still applies for a second solve
    assertEquals(s2.solveLimited(), Optional.empty());
    s2.setPropagationLimit(-1);//disable limit and solve again
    assertFalse(s2.solveLimited().get());


  }

  @Test
  public void testTimeLimit() {
    //Time limits are not implemented for mac
    Assume.assumeFalse(System.getProperty("os.name").toLowerCase().startsWith("mac"));
    Solver s2 = new Solver();
    Constraints.unsatQueens(s2,12);
    s2.setTimeLimit(1);

    long start_time = System.nanoTime();
    //s2.solve();
    assertEquals(s2.solveLimited(), Optional.empty()); //not enough time to solve
    long elapsed_seconds = (int) Math.floorDiv((System.nanoTime()-start_time),1000000000L);
    //System.out.println(elapsed_seconds);
    assertTrue(elapsed_seconds<=3);
    //limit still applies for a second solve
    assertEquals(s2.solveLimited(), Optional.empty()); //not enough time to solve
    s2.setTimeLimit(-1);//disable limit
    assertFalse(s2.solveLimited().get());
  }

  @Test
  public void testSolveIgnoresTimeLimit() {
    //Time limits are not implemented for mac
    Assume.assumeFalse(System.getProperty("os.name").toLowerCase().startsWith("mac"));
    Solver s2 = new Solver();
    Constraints.unsatQueens(s2,12);
    s2.setTimeLimit(1);
    assertFalse(s2.solve());

  }


  @Test
  public void testSecondSolverIgnoresTimeLimit() {
    //Time limits are not implemented for mac
    Assume.assumeFalse(System.getProperty("os.name").toLowerCase().startsWith("mac"));
    Solver s = new Solver();
    s.setTimeLimit(1);

    Solver s2 = new Solver();
    Constraints.unsatQueens(s2,12);
    //s2 ignores time limit
    assertFalse(s2.solveLimited().get());

  }
}
