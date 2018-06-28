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
import java.util.Optional;

import static org.junit.Assert.*;

public class ResourceLimitTest {
  @Test
  public void testConflictLimit() {

    Solver s = new Solver();

    Constraints.nqueens(s,100);
    s.setConflictLimit(-1);
    assertTrue(s.solve());

    //this requires roughly 300 conflicts to solve (with no exact guarantees)
    long n_conflicts =  s.nConflicts();
    assertTrue(n_conflicts>10);
    Solver s2 = new Solver();
    Constraints.nqueens(s2,100);
    s2.setConflictLimit(10);
    s2.setPropagationLimit(-1);
    assertEquals(s2.solveLimited(), Optional.empty());
    //limit still applies for a second solve
    assertEquals(s2.solveLimited(), Optional.empty());

    s2.setConflictLimit(-1);//disable limit and solve again
    assertEquals(s2.solveLimited().get(),true);


  }

  @Test
  public void testSolveIgnoresLimits() {
    Solver s = new Solver();
    Constraints.nqueens(s,100);
    s.setConflictLimit(10);
    s.setPropagationLimit(10);
    //solver.solve() ignores resource limitations
    assertEquals(s.solve(), true);
  }

  @Test
  public void testPropagationLimit() {

    Solver s = new Solver();

    Constraints.nqueens(s,100);
    s.setPropagationLimit(-1);
    assertTrue(s.solve());

    long n_props =  s.nPropagations();
    System.out.println(n_props);
    assertTrue(n_props>10);
    Solver s2 = new Solver();
    Constraints.nqueens(s2,100);
    s2.setPropagationLimit(10);
    s2.setConflictLimit(-1);
    assertEquals(s2.solveLimited(), Optional.empty());
    //limit still applies for a second solve
    assertEquals(s2.solveLimited(), Optional.empty());
    s2.setPropagationLimit(-1);//disable limit and solve again
    assertEquals(s2.solveLimited().get(),true);


  }

  @Test
  public void testTimeLimit() {

    Solver s = new Solver();
    Constraints.nqueens(s,100);
    long start_time = System.nanoTime();
    assertTrue(s.solve());
    long end_time = System.nanoTime();
    int elapsed_seconds = (int) Math.floorDiv((end_time-start_time),1000000000L);
    System.out.println(elapsed_seconds);
    assertTrue(elapsed_seconds>2);

    Solver s2 = new Solver();
    Constraints.nqueens(s2,100);
    s2.setTimeLimit(1);
    assertEquals(s2.solveLimited(), Optional.empty()); //not enough time to solve
    //limit still applies for a second solve
    assertEquals(s2.solveLimited(), Optional.empty()); //not enough time to solve
    s2.setTimeLimit(-1);//disable limit
    assertEquals(s2.solveLimited().get(),true);

  }

}
