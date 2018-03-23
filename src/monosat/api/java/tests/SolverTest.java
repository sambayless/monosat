import org.junit.Test;
import static org.junit.Assert.*;
import monosat.*;
import static monosat.Logic.*;


public class SolverTest {
    @Test
    public void testSolve(){
        Solver s = new Solver();
        assertEquals(s.solve(),true);
        assertEquals(s.solve(),true);
        assertEquals(s.solve(s.True()),true);
        assertEquals(s.solve(s.False()),false);
        assertEquals(s.solve(s.True(),s.False()),false);
        assertEquals(s.solve(),true);
        System.out.println("Done");
    }

    @Test
    public void testLogic(){
        Solver s = getSolver();
        Lit a = newLit();
        Lit b = newLit();
        Lit c = and(a,b);
        assertEquals(solve(),true);
    }
}
