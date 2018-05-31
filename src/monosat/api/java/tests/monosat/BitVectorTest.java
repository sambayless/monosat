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
import java.util.List;

import static org.junit.Assert.*;
/**
 * Tests the basic API functionality of bitvectors. The goal of these tests is simply to make sure
 * the JNI implementation works as expected, but does not attempt to fully test the correctness of
 * the underlying solver.
 */
public class BitVectorTest {
  static final String tricky_name =
          "~`Name-with/\"\'//<>printable_\\characters?!@#$%^&*()-+{}[]|1234567890";
  @Test
  public void getBits() {
    Solver s = new Solver();
    BitVector bv = new BitVector(s, 4);
    List<Lit> bits = bv.getBits();
    assertEquals(bits.size(), 4);
    for (Lit l : bits) {
      assertTrue(s.solve(l));
      assertTrue(s.solve(l.not()));
    }
  }

  @Test
  public void width() {
    Solver s = new Solver();
    BitVector bv = new BitVector(s, 4);
    assertEquals(bv.width(), 4);
    assertEquals(bv.size(), 4);
  }

  @Test
  public void gt() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    Lit c = bv1.gt(bv2);
    assertTrue(s.solve(c));
    long v1 = bv1.value();
    long v2 = bv2.value();
    assertTrue("BV1 > BV2", v1 > v2);
    assertTrue(s.solve(c.not()));
    v1 = bv1.value();
    v2 = bv2.value();
    assertTrue("BV1 <= BV2", v1 <= v2);
  }

  @Test
  public void geq() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    Lit c = bv1.geq(bv2);
    assertTrue(s.solve(c));
    long v1 = bv1.value();
    long v2 = bv2.value();
    assertTrue("BV1 >= BV2", v1 >= v2);
    assertTrue(s.solve(c.not()));
    v1 = bv1.value();
    v2 = bv2.value();
    assertTrue("BV1 < BV2", v1 < v2);
  }

  @Test
  public void lt() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    Lit c = bv1.lt(bv2);
    assertTrue(s.solve(c));
    long v1 = bv1.value();
    long v2 = bv2.value();
    assertTrue("BV1 < BV2", v1 < v2);
    assertTrue(s.solve(c.not()));
    v1 = bv1.value();
    v2 = bv2.value();
    assertTrue("BV1 >= BV2", v1 >= v2);
  }

  @Test
  public void leq() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    Lit c = bv1.leq(bv2);
    assertTrue(s.solve(c));
    long v1 = bv1.value();
    long v2 = bv2.value();
    assertTrue("BV1 <= BV2", v1 <= v2);
    assertTrue(s.solve(c.not()));
    v1 = bv1.value();
    v2 = bv2.value();
    assertTrue("BV1 > BV2", v1 > v2);
  }

  @Test
  public void neq() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    Lit c = bv1.neq(bv2);
    assertTrue(s.solve(c));
    long v1 = bv1.value();
    long v2 = bv2.value();
    assertTrue("BV1 != BV2", v1 != v2);
    assertTrue(s.solve(c.not()));
    v1 = bv1.value();
    v2 = bv2.value();
    assertEquals("BV1 == BV2", v1, v2);
  }

  @Test
  public void eq() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    Lit c = bv1.eq(bv2);
    assertTrue(s.solve(c));
    long v1 = bv1.value();
    long v2 = bv2.value();
    assertEquals("BV1 == BV2", v1, v2);
    assertTrue(s.solve(c.not()));
    v1 = bv1.value();
    v2 = bv2.value();
    assertTrue("BV1 != BV2", v1 != v2);
  }

  @Test
  public void eqConstAdd() {
    // Test combination of bv equality to a constant, and bv addition, which
    // once triggered a bug
    Solver s = new Solver();
    BitVector bv0 = new BitVector(s, 4);
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv3 = new BitVector(s, 4);
    monosat.Logic.assertTrue(bv3.eq(bv0.add(bv1)));
    monosat.Logic.assertTrue(bv3.eq(2));
    assertTrue(s.solve());
  }

  @Test
  public void compareBV() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);

    for (Comparison c : Comparison.values()) {
      Lit l = bv1.compare(c, bv2);
      assertTrue(s.solve(l));
      long v1 = bv1.value();
      long v2 = bv2.value();
      assertTrue(c.compare(v1, v2));
    }

    for (Comparison c : Comparison.values()) {
      Lit l = bv1.compare(c, bv2);
      assertTrue(s.solve(l.not()));
      long v1 = bv1.value();
      long v2 = bv2.value();
      assertFalse(c.compare(v1, v2));
    }
  }

  @Test
  public void compareConst() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);

    for (Comparison c : Comparison.values()) {
      // note: the extremal values (0,15) are not satisfiable for
      // some comparisons, so skipping them here
      for (long v2 = 1; v2 < 14; v2++) {
        Lit l = bv1.compare(c, v2);
        assertTrue(s.solve(l));
        long v1 = bv1.value();
        assertTrue(c.compare(v1, v2));
      }
    }
  }

  @Test
  public void unsatComparisons() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    assertFalse(s.solve(bv1.compare(Comparison.LT, 0)));
    assertFalse(s.solve(bv1.compare(Comparison.GT, 15)));
    assertTrue(s.solve());
    assertTrue(s.solve(bv1.compare(Comparison.LT, 0).not()));
    assertTrue(s.solve(bv1.compare(Comparison.GT, 15).not()));

    assertFalse(s.solve(bv1.compare(Comparison.GEQ, 0).not()));
    assertTrue(s.solve(bv1.compare(Comparison.GEQ, 0)));

    assertFalse(s.solve(bv1.compare(Comparison.LEQ, 15).not()));
    assertTrue(s.solve(bv1.compare(Comparison.LEQ, 15)));
  }

  @Test
  public void slice() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = bv1.slice(1, 3);
    assertEquals(bv2.width(), 2);
    assertFalse(s.solve(bv1.get(1), bv2.get(0).not()));
  }

  @Test
  public void concatenate() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 3);
    BitVector bv3 = bv1.concatenate(bv2);
    assertEquals(bv3.width(), 7);
    assertTrue(s.solve(bv1.get(1), bv3.get(1)));
    assertFalse(s.solve(bv1.get(1), bv3.get(1).not()));
    assertTrue(s.solve(bv2.get(1), bv3.get(5)));
    assertFalse(s.solve(bv2.get(1), bv3.get(5).not()));
  }

  @Test
  public void not() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = Logic.not(bv1);
    for (int i = 0; i < bv1.size(); i++) {
      assertFalse(s.solve(bv1.get(i), bv2.get(i)));
      assertTrue(s.solve(bv1.get(i), bv2.get(i).not()));
    }
  }

  @Test
  public void and() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.and(bv1, bv2);
    for (int i = 0; i < bv1.size(); i++) {
      assertTrue(s.solve(bv1.get(i), bv2.get(i), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i), bv2.get(i), bv3.get(i).not()));
      assertFalse(s.solve(bv1.get(i), bv2.get(i).not(), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i).not(), bv2.get(i).not(), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i).not(), bv2.get(i), bv3.get(i)));
    }
  }

  @Test
  public void nand() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.nand(bv1, bv2);
    for (int i = 0; i < bv1.size(); i++) {
      assertFalse(s.solve(bv1.get(i), bv2.get(i), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i), bv2.get(i), bv3.get(i).not()));
      assertTrue(s.solve(bv1.get(i), bv2.get(i).not(), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i).not(), bv2.get(i).not(), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i).not(), bv2.get(i), bv3.get(i)));
    }
  }

  @Test
  public void or() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.or(bv1, bv2);
    for (int i = 0; i < bv1.size(); i++) {
      assertTrue(s.solve(bv1.get(i), bv2.get(i), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i), bv2.get(i), bv3.get(i).not()));
      assertTrue(s.solve(bv1.get(i), bv2.get(i).not(), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i).not(), bv2.get(i).not(), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i).not(), bv2.get(i), bv3.get(i)));
    }
  }

  @Test
  public void nor() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.nor(bv1, bv2);
    for (int i = 0; i < bv1.size(); i++) {
      assertFalse(s.solve(bv1.get(i), bv2.get(i), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i), bv2.get(i), bv3.get(i).not()));
      assertFalse(s.solve(bv1.get(i), bv2.get(i).not(), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i).not(), bv2.get(i).not(), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i).not(), bv2.get(i), bv3.get(i)));
    }
  }

  @Test
  public void xor() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.xor(bv1, bv2);
    for (int i = 0; i < bv1.size(); i++) {
      assertFalse(s.solve(bv1.get(i), bv2.get(i), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i), bv2.get(i), bv3.get(i).not()));
      assertTrue(s.solve(bv1.get(i), bv2.get(i).not(), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i).not(), bv2.get(i).not(), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i).not(), bv2.get(i), bv3.get(i)));
    }
  }

  @Test
  public void xnor() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.xnor(bv1, bv2);
    for (int i = 0; i < bv1.size(); i++) {
      assertTrue(s.solve(bv1.get(i), bv2.get(i), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i), bv2.get(i), bv3.get(i).not()));
      assertFalse(s.solve(bv1.get(i), bv2.get(i).not(), bv3.get(i)));
      assertTrue(s.solve(bv1.get(i).not(), bv2.get(i).not(), bv3.get(i)));
      assertFalse(s.solve(bv1.get(i).not(), bv2.get(i), bv3.get(i)));
    }
  }

  @Test
  public void add() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.add(bv1, bv2);
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        assertTrue(s.solve(bv1.eq(i), bv2.eq(j)));
        assertEquals(bv1.value(), i);
        assertEquals(bv2.value(), j);
        assertEquals(bv3.value(), i + j);
        assertFalse(s.solve(bv1.eq(i), bv2.eq(j), bv3.neq(i + j)));
      }
    }
  }

  @Test
  public void subtract() {
    Solver s = new Solver();
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = bv1.subtract(bv2);
    for (int i = 0; i < 7; i++) {
      for (int j = 0; j <= i; j++) {
        assertTrue(s.solve(bv1.eq(i), bv2.eq(j)));
        assertEquals(bv1.value(), i);
        assertEquals(bv2.value(), j);
        assertEquals(bv3.value(), i - j);
        assertFalse(s.solve(bv1.eq(i), bv2.eq(j), bv3.neq(i - j)));
      }
    }
  }

  @Test
  public void getBitVectors() {
    Solver s = new Solver();
    assertEquals(s.getBitVectors().size(), 0);
    BitVector bv1 = new BitVector(s, 4);
    BitVector bv2 = new BitVector(s, 4);
    BitVector bv3 = Logic.add(bv1, bv2);
    assertEquals(s.getBitVectors().size(), 3);
    s.assertTrue(bv1.gt(1));
    s.assertTrue(bv2.gt(bv1));
    s.solve();
    for (BitVector bv : s.getBitVectors()) {
      long val = bv.value();
      assert (val > 0);
    }
    assertEquals(s.getBitVectors().size(), 3);
  }

  @Test
  public void maxMinBVSize() {
    Solver s = new Solver();
    BitVector bv0 = new BitVector(s, 1);
    assertEquals(bv0.width(),1);

    BitVector bv1 = new BitVector(s, 64);
    assertEquals(bv1.width(),64);

    try {
      new BitVector(s, 65,"MyBitvector");
      fail("Maximum BV size is 64");
    } catch (IllegalArgumentException except) {
      // ok
    }
    try {
      new BitVector(s, 0,"MyBitvector");
      fail("Minimum BV size is 0");
    } catch (IllegalArgumentException except) {
      // ok
    }

  }


  @Test
  public void testNamedBVs() {
    monosat.Solver s = new monosat.Solver();
    BitVector bv0 = new BitVector(s, 4);
    BitVector bv1 = new BitVector(s, 4,"");
    BitVector bv2 = new BitVector(s, 4,"");

    BitVector bv3 = new BitVector(s, 4,"MyBitvector");

    try {
      new BitVector(s, 4,"MyBitvector");
      fail("No two bvs can have the same name");
    } catch (IllegalArgumentException except) {
      // ok
    }

    try {
      new BitVector(s, 4,"Name With Spaces");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }
    BitVector bv4 =  new BitVector(s, 4,tricky_name);


    try {
      new BitVector(s, 4,"Name With \n NewLine");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }
    try {
      new BitVector(s, 4,"Name With \t tab");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }

    assertEquals(bv0.name(), "");
    assertEquals(bv1.name(), "");
    assertEquals(bv2.name(), "");
    assertEquals(bv3.name(), "MyBitvector");
    assertEquals(bv4.name(), tricky_name);
    try {
      s.getBitvector("");
      fail("Expected a bad name exception");
    } catch (IllegalArgumentException except) {
      // ok
    }

    assertEquals(s.getBitvector("MyBitvector"), bv3);
    assertEquals(s.getBitvector(tricky_name), bv4);
  }
      @Test
  public void testLoadingBV() throws IOException {
    File file = File.createTempFile("test", ".gnf");
    String filename = file.getAbsolutePath().toString();
    file.delete();

    {
      monosat.Solver s = new monosat.Solver();
      BitVector bv0 = new BitVector(s, 4);
      BitVector bv1 = new BitVector(s, 4,"");
      BitVector bv2 = new BitVector(s, 4,"");

      BitVector bv3 = new BitVector(s, 4,"MyBitvector");

      try {
        new BitVector(s, 4,"MyBitvector");
        fail("No two bvs can have the same name");
      } catch (IllegalArgumentException except) {
        // ok
      }

      try {
        new BitVector(s, 4,"Name With Spaces");
        fail("Expected a bad name exception");
      } catch (IllegalArgumentException except) {
        // ok
      }
      BitVector bv4 =  new BitVector(s, 4,tricky_name);


      try {
        new BitVector(s, 4,"Name With \n NewLine");
        fail("Expected a bad name exception");
      } catch (IllegalArgumentException except) {
        // ok
      }
      try {
        new BitVector(s, 4,"Name With \t tab");
        fail("Expected a bad name exception");
      } catch (IllegalArgumentException except) {
        // ok
      }

      assertEquals(bv0.name(), "");
      assertEquals(bv1.name(), "");
      assertEquals(bv2.name(), "");
      assertEquals(bv3.name(), "MyBitvector");
      assertEquals(bv4.name(), tricky_name);
      try {
        s.getBitvector("");
        fail("Expected a bad name exception");
      } catch (IllegalArgumentException except) {
        // ok
      }

      assertEquals(s.getBitvector("MyBitvector"), bv3);
      assertEquals(s.getBitvector(tricky_name), bv4);

      s.close();
    }

    monosat.Solver s = new monosat.Solver();
    assertTrue(s.solve());
    s.loadConstraints(filename);
    assertTrue(s.solve());
    BitVector bv1 = s.getBitvector("MyBitvector");
    assertEquals(bv1.width(),4);
    BitVector bv2 = s.getBitvector(tricky_name);

    assertEquals(s.getBitvector("MyBitvector"), bv1);
    assertEquals(s.getBitvector(tricky_name), bv2);

    assertTrue(s.solve(bv1.eq( 1)));
  }
}
