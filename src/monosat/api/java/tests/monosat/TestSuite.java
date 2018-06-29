package monosat;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

@RunWith(Suite.class)
@Suite.SuiteClasses({
  monosat.BitVectorTest.class,
  monosat.GraphTest.class,
  monosat.LogicTest.class,
  monosat.SolverTest.class,
  monosat.LitTest.class,
  monosat.PseudoBooleanTest.class,
  monosat.ResourceLimitTest.class
})
public class TestSuite {}
