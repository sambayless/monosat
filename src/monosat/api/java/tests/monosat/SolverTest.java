package monosat;

import org.junit.Test;
import static org.junit.Assert.*;
import monosat.*;


public class SolverTest {
    @Test
    public void testSolve(){
        Solver s = new Solver();
        assertEquals(s.solve(),true);
        assertEquals(s.solve(),true);
        assertEquals(s.solve(s.getTrue()),true);
        assertEquals(s.solve(s.getFalse()),false);
        assertEquals(s.solve(s.getTrue(),s.getFalse()),false);
        assertEquals(s.solve(),true);
        System.out.println("Done");
    }

}
