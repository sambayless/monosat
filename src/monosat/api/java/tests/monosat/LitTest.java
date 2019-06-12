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
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.*;

public class LitTest {
  static final String tricky_name =
      "`Name-with/\"\'//<>printable_\\characters?!@#$%^&*()-+{}[]|1234567890";

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

    assertEquals(Lit.True.toVar(),Lit.False.toVar());
    assertEquals(Lit.True.toVar(),0);
    assertEquals(Lit.True.toVar()+1,a.toVar());
    assertEquals(Lit.True.toVar()+2,b.toVar());
    assertEquals(Lit.True.toVar()+3,c.toVar());

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

      assertEquals(a.name(), a.not().name());
      assertEquals(b.name(), b.not().name());
      assertEquals("~" + c.name(),c.not().name());
      assertEquals("~" +e.name(), e.not().name());

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
    public void testRenamingLits() {
        monosat.Solver s = new monosat.Solver();
        monosat.Lit a = new monosat.Lit(s);
        monosat.Lit b = new monosat.Lit(s, "");

        monosat.Lit c = new monosat.Lit(s, "MyLiteral");


        s.addName(c,"MyLiteral1");

        assertEquals(a.name(),"");
        try{
            s.addName(a,"MyLiteral1");
            fail("No two variables can have the same name");
        } catch (IllegalArgumentException except) {
            // ok
        }


        s.addName(a,"MyLiteral2");
        assertEquals(a.name(),"MyLiteral2");
        s.addName(a,"");
        assertEquals(a.name(),"MyLiteral2");
        s.addName(a,"MyLiteral4");
        //After adding a second name, the first name of the literal does not change
        assertEquals(a.name(),"MyLiteral2");

        try {
            s.addName(b,"MyLiteral");
            fail("No two variables can have the same name");
        } catch (IllegalArgumentException except) {
            // ok
        }

        assertEquals(b.name(),"");


        try {
            s.addName(b,"MyLiteral2");
            fail("No two variables can have the same name");
        } catch (IllegalArgumentException except) {
            // ok
        }

        s.addName(b.not(),"MyLiteral3");


        assertEquals(a.name(), "MyLiteral2");
        assertEquals(b.name(), "~MyLiteral3");
        assertEquals(c.name(), "MyLiteral");


        assertEquals("~"+a.name(), a.not().name());
        assertEquals(b.name(),"~"+ b.not().name());
        assertEquals("~"+c.name(),c.not().name());


        try {
            s.getLiteral("");
            fail("Expected a bad name exception");
        } catch (IllegalArgumentException except) {
            // ok
        }

        assertEquals(s.getLiteral("True"), Lit.True); // "True" is always named in the solver
        assertEquals(s.getLiteral("False"), Lit.False); // "False" is always named in the solver
        assertEquals(s.getLiteral("MyLiteral"), c);
        assertEquals(s.getLiteral("MyLiteral2"), a);
        assertEquals(s.getLiteral("MyLiteral3"), b.not());

        assertTrue(a.hasName("MyLiteral2"));
        assertFalse(a.hasName("~MyLiteral2"));
        assertFalse(a.not().hasName("MyLiteral2"));
        assertTrue(a.not().hasName("~MyLiteral2"));
        assertTrue(c.hasName("MyLiteral"));
        assertFalse(c.hasName("~MyLiteral"));
        assertTrue(c.not().hasName("~MyLiteral"));
        assertTrue(b.hasName("~MyLiteral3"));
        assertTrue(b.not().hasName("MyLiteral3"));
        assertFalse(a.hasName("MyLiteral3"));
        assertFalse(a.hasName("MyLiteral"));
        assertFalse(c.hasName("MyLiteral2"));
        assertFalse(b.hasName("MyLiteral"));
        assertFalse(b.hasName("MyLiteral2"));

        assertTrue(a.hasName("MyLiteral4"));
        assertTrue(a.not().hasName("~MyLiteral4"));

        ArrayList<String> anames = new ArrayList<>();
        for(String str:a.names()){
            anames.add(str);
            assertTrue(a.hasName(str));
            assertFalse(b.hasName(str));
            assertFalse(c.hasName(str));
        }
        assertTrue(anames.size()==2);
        assertTrue(anames.contains("MyLiteral2"));
        assertTrue(anames.contains("MyLiteral4"));

        for(String str:b.names()){
            assertTrue(b.hasName(str));
            assertFalse(a.hasName(str));
            assertFalse(c.hasName(str));
        }
        for(String str:c.names()){
            assertTrue(c.hasName(str));
            assertFalse(a.hasName(str));
            assertFalse(b.hasName(str));
        }

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

      s.closeConstraintFile();
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
    Lit n2 = new Lit(s2, "MyLiteral2");

    s2.loadConstraints(filename);
    assertTrue(s2.solve());
    Lit c2 = s2.getLiteral("MyLiteral");
    Lit e2 = s2.getLiteral(tricky_name);
    Lit m2 = new Lit(s2, "MyLiteral3");

    assertEquals(s2.getLiteral("True"), Lit.True); // "True" is always named in the solver
    assertEquals(s2.getLiteral("False"), Lit.False); // "False" is always named in the solver
    assertEquals(s2.getLiteral("MyLiteral"), c2);
    assertEquals(s2.getLiteral(tricky_name), e2);

    assertTrue(s2.solve(c2));
    assertTrue(s2.solve(e2));
    assertFalse(s2.solve(c2, e2));

    // The solver should (now) maintain integer mappings of literals after loading from disk
    assertEquals(c2.toInt(), c.toInt());
    assertEquals(e2.toInt(), e.toInt());
  }

    @Test
    public void testLoadingNegatedLits() throws IOException {
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
        Lit n2 = new Lit(s2, "MyLiteral2");

        s2.loadConstraints(filename);
        assertTrue(s2.solve());
        Lit c2 = s2.getLiteral("MyLiteral");
        Lit c2neg = s2.getLiteral("~MyLiteral");
        Lit e2 = s2.getLiteral(tricky_name);
        Lit e2neg = s2.getLiteral("~"+tricky_name);
        Lit m2neg = new Lit(s2, "MyLiteral3");

        assertEquals(s2.getLiteral("True"), Lit.True); // "True" is always named in the solver
        assertEquals(s2.getLiteral("False"), Lit.False); // "False" is always named in the solver
        assertEquals(s2.getLiteral("MyLiteral"), c2);
        assertEquals(s2.getLiteral(tricky_name), e2);

        assertEquals(c2.not(),c2neg);
        assertEquals(e2.not(),e2neg);

    }

    @Test
    public void testLitIterator() {
        monosat.Solver s = new monosat.Solver();
        assertEquals(1,s.nVars());
        monosat.Lit a = new monosat.Lit(s);
        assertEquals(2,s.nVars());
        monosat.Lit b = new monosat.Lit(s, "");
        assertEquals(3,s.nVars());
        monosat.Lit b2 = new monosat.Lit(s, ""); // it is ok for multiple literals to have empty names
        monosat.Lit c = new monosat.Lit(s, "MyLiteral");
        assertEquals(5,s.nVars());
        try {
            monosat.Lit c2 = new monosat.Lit(s, "MyLiteral");
            fail("No two variables can have the same name");
        } catch (IllegalArgumentException except) {
            // ok
        }
        assertEquals(5,s.nVars());
        try {
            monosat.Lit t = new monosat.Lit(s, "True");
            fail("No two variables can have the same name");
        } catch (IllegalArgumentException except) {
            // ok
        }
        assertEquals(5,s.nVars());
        try {
            monosat.Lit t = new monosat.Lit(s, "False");
            fail("No two variables can have the same name");
        } catch (IllegalArgumentException except) {
            // ok
        }
        assertEquals(5,s.nVars());
        try {
            monosat.Lit d = new monosat.Lit(s, "Name With Spaces");
            fail("Expected a bad name exception");
        } catch (IllegalArgumentException except) {
            // ok
        }

        monosat.Lit e = new monosat.Lit(s, tricky_name);
        assertEquals(6,s.nVars());
        try {
            monosat.Lit f = new monosat.Lit(s, "Name With \n NewLine");
            fail("Expected a bad name exception");
        } catch (IllegalArgumentException except) {
            // ok
        }
        assertEquals(6,s.nVars());
        try {
            monosat.Lit g = new monosat.Lit(s, "Name With \t tab");
            fail("Expected a bad name exception");
        } catch (IllegalArgumentException except) {
            // ok
        }
        assertEquals(6,s.nVars());
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

        assertEquals(6,s.nVars());
        {
          Iterator<Lit> it = s.literals().iterator();
          assertEquals(it.next(), Lit.True);
          assertEquals(it.next(), a);
          assertEquals(it.next(), b);
          assertEquals(it.next(), b2);
          assertEquals(it.next(), c);
          assertEquals(it.next(), e);
          assertFalse(it.hasNext());
          try {
            it.next();
            fail("Expected out of bounds exception");
          } catch (IndexOutOfBoundsException except) {
            // ok
          }
       }
        {
          Iterator<Lit> it2 = s.namedLiterals().iterator();
          assertEquals(it2.next(), Lit.True);
          assertEquals(it2.next(), c);
          assertEquals(it2.next(), e);
          assertFalse(it2.hasNext());
          try {
            it2.next();
            fail("Expected out of bounds exception");
          } catch (IndexOutOfBoundsException except) {
            // ok
          }
        }
    }


    @Test
    public void testLoadingLitsIterators() throws IOException {
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
            assertEquals(6,s.nVars());
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
            assertEquals(6,s.nVars());
            assertTrue(s.solve());
        }

        monosat.Solver s2 = new monosat.Solver();
        assertTrue(s2.solve());
        assertEquals(1,s2.nVars());
        Lit n2 = new Lit(s2, "MyLiteral2");
        assertEquals(2,s2.nVars());
        s2.loadConstraints(filename); //Note that literal n2, the first declared variable (after 'True'), will
        //be mapped to the same literal as literal 'a' in the GNF formula.
        assertEquals(6,s2.nVars());
        assertTrue(s2.solve());
        Lit c2 = s2.getLiteral("MyLiteral");
        Lit e2 = s2.getLiteral(tricky_name);
        Lit m2 = new Lit(s2, "MyLiteral3");
        assertEquals(7,s2.nVars());
        assertEquals(s2.getLiteral("True"), Lit.True); // "True" is always named in the solver
        assertEquals(s2.getLiteral("False"), Lit.False); // "False" is always named in the solver
        assertEquals(s2.getLiteral("MyLiteral"), c2);
        assertEquals(s2.getLiteral(tricky_name), e2);

        assertTrue(s2.solve(c2));
        assertTrue(s2.solve(e2));
        assertFalse(s2.solve(c2, e2));

        // The solver should (now) maintain integer mappings of literals after loading from disk

        {
            Iterator<Lit> it = s2.literals().iterator();
            assertEquals(it.next(), Lit.True);
            assertEquals(it.next(), n2);
            it.next();
            it.next();
            assertEquals(it.next(), c2);
            assertEquals(it.next(), e2);
            assertEquals(it.next(), m2);
            assertFalse(it.hasNext());
            try {
                it.next();
                fail("Expected out of bounds exception");
            } catch (IndexOutOfBoundsException except) {
                // ok
            }

        }
        {
            Iterator<Lit> it2 = s2.namedLiterals().iterator();
            assertEquals(it2.next(), Lit.True);
            assertEquals(it2.next(), n2);
            assertEquals(it2.next(), c2);
            assertEquals(it2.next(), e2);
            assertEquals(it2.next(), m2);
            assertFalse(it2.hasNext());
            try {
                it2.next();
                fail("Expected out of bounds exception");
            } catch (IndexOutOfBoundsException except) {
                // ok
            }
        }
    }


    @Test
    public void testLoadingLitsNameClash() throws IOException {
        File file = File.createTempFile("test", ".gnf");
        String filename = file.getAbsolutePath().toString();
        file.delete();
        {
            monosat.Solver s = new monosat.Solver("", filename);
            monosat.Lit a = new monosat.Lit(s,"MyLiteral1");
            s.flushConstraintFile();
        }

        monosat.Solver s2 = new monosat.Solver();
        assertTrue(s2.solve());
        assertEquals(1,s2.nVars());
        Lit a2 = new Lit(s2, "MyLiteral2");

        //If you load constraints, and those constraints rename an existing, already
        //named variable, then that is no longer an exception, as multiple names
        //per literal are now supported.
        s2.loadConstraints(filename);

        assertTrue(a2.hasName("MyLiteral2"));
        assertTrue(a2.hasName("MyLiteral1"));

    }


    @Test
    public void testComparable() {
        try(monosat.Solver s = new monosat.Solver()) {
            List<Lit> lits = new ArrayList<>();
            for(int i = 0;i<10;i++){
                lits.add(new Lit(s));
            }
            Collections.shuffle(lits);
            Collections.sort(lits);
        }
    }
}
