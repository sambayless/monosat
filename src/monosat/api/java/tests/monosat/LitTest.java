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

import static org.junit.Assert.*;

public class LitTest {
  static final String tricky_name =
      "~`Name-with/\"\'//<>printable_\\characters?!@#$%^&*()-+{}[]|1234567890";

  @Test
  public void testGlobals() {
    monosat.Solver s = new monosat.Solver();
    assertTrue(s.solve());
    assertTrue(s.solve());
    assertTrue(s.solve(monosat.Lit.True));
    assertFalse(s.solve(monosat.Lit.False));
    assertFalse(s.solve(monosat.Lit.True, monosat.Lit.False));
    assertTrue(s.solve());
  }

  @Test
  public void testNewLits() {
    monosat.Solver s = new monosat.Solver();
    monosat.Lit a = new monosat.Lit(s);
    monosat.Lit b = new monosat.Lit(s);
    monosat.Lit c = new monosat.Lit(s);

    assertTrue(s.solve(a));
    assertTrue(a.value());
    assertTrue(s.solve(b.not()));
    assertFalse(b.value());
    assertTrue(s.solve(a, b.not()));
    assertTrue(a.value());
    assertFalse(b.value());
  }

  @Test
  public void testMaybeValues() {
    monosat.Solver s = new monosat.Solver();
    monosat.Lit a = new monosat.Lit(s);
    monosat.Lit b = new monosat.Lit(s);
    monosat.Lit c = new monosat.Lit(s);

    assertFalse(a.possibleValue().isPresent());

    assertTrue(s.solve(a));
    assertTrue(a.possibleValue().isPresent());
    assertTrue(s.solve(b.not()));
    assertTrue(a.possibleValue().isPresent());
    assertTrue(b.possibleValue().isPresent());
    assertTrue(s.solve(a, b.not()));
    assertTrue(a.possibleValue().isPresent());
    assertTrue(b.possibleValue().isPresent());
  }

  @Test
  public void testDefaultValues() {
    monosat.Solver s = new monosat.Solver();
    monosat.Lit a = new monosat.Lit(s);

    assertFalse(a.value(false));
    assertTrue(a.value(true));

    assertTrue(s.solve(a));
    assertTrue(a.value(false));
    assertTrue(a.value(true));
  }

  @Test
  public void testConstLits() {
    monosat.Solver s = new monosat.Solver();
    monosat.Lit a = new monosat.Lit(s);
    monosat.Lit b = new monosat.Lit(s);
    monosat.Lit c = new monosat.Lit(s);
    assertTrue(monosat.Lit.True.isConst());
    assertTrue(monosat.Lit.True.isConstTrue());
    assertFalse(monosat.Lit.True.isConstFalse());

    assertFalse(a.isConst());
    s.addClause(a);
    assertTrue(a.isConst());
    assertTrue(a.isConstTrue());
    assertFalse(a.isConstFalse());

    s.addClause(b.not());
    assertTrue(b.isConst());
    assertFalse(b.isConstTrue());
    assertTrue(b.isConstFalse());
  }

  @Test
  public void testNamedLits() {
    monosat.Solver s = new monosat.Solver();
    monosat.Lit a = new monosat.Lit(s);
    monosat.Lit b = new monosat.Lit(s, "");
    monosat.Lit b2 = new monosat.Lit(s, ""); // it is ok for multiple literals to have empty names
    monosat.Lit c = new monosat.Lit(s, "MyLiteral");

    try {
      monosat.Lit c2 = new monosat.Lit(s, "MyLiteral");
      fail("No two variables can have the same name");
    } catch (IllegalArgumentException except) {
      // ok
    }

    try {
      monosat.Lit t = new monosat.Lit(s, "True");
      fail("No two variables can have the same name");
    } catch (IllegalArgumentException except) {
      // ok
    }

    try {
      monosat.Lit t = new monosat.Lit(s, "False");
      fail("No two variables can have the same name");
    } catch (IllegalArgumentException except) {
      // ok
    }

    try {
      monosat.Lit d = new monosat.Lit(s, "Name With Spaces");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }

    monosat.Lit e = new monosat.Lit(s, tricky_name);

    try {
      monosat.Lit f = new monosat.Lit(s, "Name With \n NewLine");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }
    try {
      monosat.Lit g = new monosat.Lit(s, "Name With \t tab");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }

    assertEquals(a.name(), "");
    assertEquals(b.name(), "");
    assertEquals(c.name(), "MyLiteral");
    assertEquals(e.name(), tricky_name);

    try {
      s.getLiteral("");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }

    assertEquals(s.getLiteral("True"), Lit.True); // "True" is always named in the solver
    assertEquals(s.getLiteral("False"), Lit.False); // "False" is always named in the solver
    assertEquals(s.getLiteral("MyLiteral"), c);
    assertEquals(s.getLiteral(tricky_name), e);
  }

  @Test
  public void testLoadingLits() throws IOException {
    File file = File.createTempFile("test", ".gnf");
    String filename = file.getAbsolutePath().toString();
    file.delete();

    {
      monosat.Solver s = new monosat.Solver("", filename);
      monosat.Lit a = new monosat.Lit(s);
      monosat.Lit b = new monosat.Lit(s, "");
      monosat.Lit b2 = new monosat.Lit(s, ""); // it is ok for multiple literals to have empty names
      monosat.Lit c = new monosat.Lit(s, "MyLiteral");

      monosat.Lit e = new monosat.Lit(s, tricky_name);

      s.addClause(c.not(), e.not());

      try {
        monosat.Lit f = new monosat.Lit(s, "Name With \n NewLine");
        fail("Expected a bad name exception");
      } catch (IllegalArgumentException except) {
        // ok
      }
      try {
        monosat.Lit g = new monosat.Lit(s, "Name With \t tab");
        fail("Expected a bad name exception");
      } catch (IllegalArgumentException except) {
        // ok
      }

      s.close();
    }

    monosat.Solver s = new monosat.Solver();
    assertTrue(s.solve());
    s.loadConstraints(filename);
    assertTrue(s.solve());
    Lit c = s.getLiteral("MyLiteral");
    Lit e = s.getLiteral(tricky_name);

    try {
      monosat.Lit f = new monosat.Lit(s, "MyLiteral");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }

    assertEquals(s.getLiteral("True"), Lit.True); // "True" is always named in the solver
    assertEquals(s.getLiteral("False"), Lit.False); // "False" is always named in the solver
    assertEquals(s.getLiteral("MyLiteral"), c);
    assertEquals(s.getLiteral(tricky_name), e);

    assertTrue(s.solve(c));
    assertTrue(s.solve(e));
    assertFalse(s.solve(c, e));
  }


  @Test
  public void testLoadingLitsIntegers() throws IOException {
    File file = File.createTempFile("test", ".gnf");
    String filename = file.getAbsolutePath().toString();
    file.delete();


    monosat.Solver s = new monosat.Solver("", filename);
    monosat.Lit a = new monosat.Lit(s);
    monosat.Lit b = new monosat.Lit(s, "");
    monosat.Lit b2 = new monosat.Lit(s, ""); // it is ok for multiple literals to have empty names
    monosat.Lit c = new monosat.Lit(s, "MyLiteral");

    monosat.Lit e = new monosat.Lit(s, tricky_name);

    s.addClause(c.not(), e.not());

    try {
      monosat.Lit f = new monosat.Lit(s, "Name With \n NewLine");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }
    try {
      monosat.Lit g = new monosat.Lit(s, "Name With \t tab");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }

    assertTrue(s.solve());


    monosat.Solver s2 = new monosat.Solver();
    assertTrue(s2.solve());
    Lit n2 = new Lit(s2,"MyLiteral2");

    s2.loadConstraints(filename);
    assertTrue(s2.solve());
    Lit c2 = s2.getLiteral("MyLiteral");
    Lit e2 = s2.getLiteral(tricky_name);
    Lit m2 = new Lit(s2,"MyLiteral3");

    assertEquals(s2.getLiteral("True"), Lit.True); // "True" is always named in the solver
    assertEquals(s2.getLiteral("False"), Lit.False); // "False" is always named in the solver
    assertEquals(s2.getLiteral("MyLiteral"), c2);
    assertEquals(s2.getLiteral(tricky_name), e2);

    assertTrue(s2.solve(c2));
    assertTrue(s2.solve(e2));
    assertFalse(s2.solve(c2, e2));

    //The solver should (now) maintain integer mappings of literals after loading from disk
    assertEquals(c2.toInt(),c.toInt());
    assertEquals(e2.toInt(),e.toInt());
  }
}
