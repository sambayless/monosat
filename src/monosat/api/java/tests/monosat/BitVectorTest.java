/***************************************************************************************[Solver.cc]
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
 **************************************************************************************************/

package monosat;

import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
/**
 * Tests the basic API functionality of bitvectors.
 * The goal of these tests is simply to make sure the JNI implementation works as expected,
 * but does not attempt to fully test the correctness of the underlying solver.
 */
public class BitVectorTest {

    @Test
    public void getBits() {
        Solver s = new Solver();
        BitVector bv = new BitVector(s,4);
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
        BitVector bv = new BitVector(s,4);
        assertEquals(bv.width(), 4);
        assertEquals(bv.size(), 4);
    }

    @Test
    public void gt() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        Lit c = bv1.gt(bv2);
        assertTrue(s.solve(c));
        long v1 = s.getValue(bv1);
        long v2 = s.getValue(bv2);
        assertTrue("BV1 > BV2", v1 > v2);
        assertTrue(s.solve(c.not()));
        v1 = s.getValue(bv1);
        v2 = s.getValue(bv2);
        assertTrue("BV1 <= BV2", v1 <= v2);
    }

    @Test
    public void geq() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        Lit c = bv1.geq(bv2);
        assertTrue(s.solve(c));
        long v1 = s.getValue(bv1);
        long v2 = s.getValue(bv2);
        assertTrue("BV1 >= BV2", v1 >= v2);
        assertTrue(s.solve(c.not()));
        v1 = s.getValue(bv1);
        v2 = s.getValue(bv2);
        assertTrue("BV1 < BV2", v1 < v2);
    }

    @Test
    public void lt() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        Lit c = bv1.lt(bv2);
        assertTrue(s.solve(c));
        long v1 = s.getValue(bv1);
        long v2 = s.getValue(bv2);
        assertTrue("BV1 < BV2", v1 < v2);
        assertTrue(s.solve(c.not()));
        v1 = s.getValue(bv1);
        v2 = s.getValue(bv2);
        assertTrue("BV1 >= BV2", v1 >= v2);
    }

    @Test
    public void leq() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        Lit c = bv1.leq(bv2);
        assertTrue(s.solve(c));
        long v1 = s.getValue(bv1);
        long v2 = s.getValue(bv2);
        assertTrue("BV1 <= BV2", v1 <= v2);
        assertTrue(s.solve(c.not()));
        v1 = s.getValue(bv1);
        v2 = s.getValue(bv2);
        assertTrue("BV1 > BV2", v1 > v2);
    }

    @Test
    public void neq() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        Lit c = bv1.neq(bv2);
        assertTrue(s.solve(c));
        long v1 = s.getValue(bv1);
        long v2 = s.getValue(bv2);
        assertTrue("BV1 != BV2", v1 != v2);
        assertTrue(s.solve(c.not()));
        v1 = s.getValue(bv1);
        v2 = s.getValue(bv2);
        assertTrue("BV1 == BV2", v1 == v2);
    }

    @Test
    public void eq() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        Lit c = bv1.eq(bv2);
        assertTrue(s.solve(c));
        long v1 = s.getValue(bv1);
        long v2 = s.getValue(bv2);
        assertTrue("BV1 == BV2", v1 == v2);
        assertTrue(s.solve(c.not()));
        v1 = s.getValue(bv1);
        v2 = s.getValue(bv2);
        assertTrue("BV1 != BV2", v1 != v2);
    }


    @Test
    public void eqConstAdd() {
        //Test combination of bv equality to a constant, and bv addition, which
        //once triggered a bug
        Solver s = new Solver();
        BitVector bv0 = new BitVector(s,4);
        BitVector bv1 = new BitVector(s,4);
        BitVector bv3 = new BitVector(s,4);
        monosat.Logic.assertTrue(bv3.eq (bv0.add(bv1)));
        monosat.Logic.assertTrue(bv3.eq (2));
        assertTrue(s.solve());
    }



    @Test
    public void slice() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = bv1.slice(1, 3);
        assertEquals(bv2.width(), 2);
        assertFalse(s.solve(bv1.get(1), bv2.get(0).not()));
    }

    @Test
    public void append() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,3);
        BitVector bv3 = bv1.append(bv2);
        assertEquals(bv3.width(), 7);
        assertTrue(s.solve(bv1.get(1), bv3.get(1)));
        assertFalse(s.solve(bv1.get(1), bv3.get(1).not()));
        assertTrue(s.solve(bv2.get(1), bv3.get(5)));
        assertFalse(s.solve(bv2.get(1), bv3.get(5).not()));
    }

    @Test
    public void not() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = Logic.not(bv1);
        for (int i = 0; i < bv1.size(); i++) {
            assertFalse(s.solve(bv1.get(i), bv2.get(i)));
            assertTrue(s.solve(bv1.get(i), bv2.get(i).not()));
        }
    }

    @Test
    public void and() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = Logic.and(bv1,bv2);
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
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = Logic.nand(bv1,bv2);
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
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = Logic.or(bv1,bv2);
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
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = Logic.nor(bv1,bv2);
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
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = Logic.xor(bv1,bv2);
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
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = Logic.xnor(bv1,bv2);
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
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = Logic.add(bv1,bv2);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertTrue(s.solve(bv1.eq(i), bv2.eq(j)));
                assertEquals(s.getValue(bv1), i);
                assertEquals(s.getValue(bv2), j);
                assertEquals(s.getValue(bv3), i + j);
                assertFalse(s.solve(bv1.eq(i), bv2.eq(j), bv3.neq(i + j)));
            }
        }
    }

    @Test
    public void subtract() {
        Solver s = new Solver();
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);
        BitVector bv3 = bv1.subtract(bv2);
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j <= i; j++) {
                assertTrue(s.solve(bv1.eq(i), bv2.eq(j)));
                assertEquals(s.getValue(bv1), i);
                assertEquals(s.getValue(bv2), j);
                assertEquals(s.getValue(bv3), i - j);
                assertFalse(s.solve(bv1.eq(i), bv2.eq(j), bv3.neq(i - j)));
            }
        }
    }
}