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

import java.util.ArrayList;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class GraphTest {

    @Test
    public void nNodes() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 10; i++) {
            g.addNode();
            assertEquals(g.nNodes(), i + 1);
        }
    }

    @Test
    public void nEdges() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
        /*

        0 *--* 1
          |/ |
        2 *--* 3

         */
        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2);
        Lit e_2_3 = g.addEdge(2, 3);

        assertEquals(g.nEdges(), 5);
    }

    @Test
    public void addNode() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 10; i++) {
            g.addNode();
            assertEquals(g.nNodes(), i + 1);
        }
    }

    @Test
    public void bitwidth() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        assertEquals(g.bitwidth(), -1);
        Graph g2 = new Graph(s, 4);
        assertEquals(g2.bitwidth(), 4);
    }

    @Test
    public void addEdge() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
        /*

        0 *--* 1
          |/ |
        2 *--* 3

         */
        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2);
        Lit e_2_3 = g.addEdge(2, 3);

        assertEquals(s.solve(e_0_1, e_0_2.negate(), e_1_3), true);
    }

    @Test
    public void reaches() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
        /*

        0 *--* 1
          |/ |
        2 *--* 3

         */
        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2);
        Lit e_2_3 = g.addEdge(2, 3);

        Lit r = g.reaches(0, 3);
        assertTrue(s.solve(r));
        assertFalse(s.solve(r, e_0_1.negate(), e_2_3.negate()));
        assertTrue(s.solve(r, e_0_2.negate(), e_1_3.negate()));

        assertFalse(s.solve(r, e_0_2.negate(), e_1_3.negate(), e_2_3.negate()));

        assertTrue(s.solve(r, e_0_2.negate(), e_1_3.negate()));
        //There should only be one solution to this: 0->1, 1->2, 2->3
        ArrayList<Integer> nodes = g.getPathNodes(r);
        ArrayList<Lit> edges = g.getPathEdges(r);
        assertEquals(edges.size(), 3);
        assertEquals(nodes.size(), 4);

        assertEquals(edges.get(0), e_0_1);
        assertEquals(edges.get(1), e_1_2);
        assertEquals(edges.get(2), e_2_3);
        assertTrue(s.getValue(e_0_1));
        assertTrue(s.getValue(e_1_2));
        assertTrue(s.getValue(e_2_3));
        assertFalse(s.getValue(e_0_2));
        assertFalse(s.getValue(e_1_3));

        assertTrue(s.solve(r));

    }


    @Test
    public void maximumFlow_geq() {
        Solver s = new Solver();
        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*
  
          0 *--* 1
            |/ |
          2 *--* 3
  
           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector cmp = s.bv(4);
        Lit f = g.compareMaximumFlow(0, 3, cmp, Comparison.GEQ);
        assertTrue(s.solve(f));
        assertFalse(s.solve(f, cmp.eq(3)));
        assertTrue(s.solve(f, cmp.eq(2)));
        assertTrue(s.solve(f, cmp.eq(1)));
        assertTrue(s.solve(f, e_0_1.negate(), e_2_3.negate()));
        assertFalse(s.solve(f, e_0_1.negate(), e_2_3.negate(), cmp.eq(1)));
        assertTrue(s.solve(f, e_0_2.negate(), e_1_3.negate()));
        assertTrue(s.solve(f));
    }

    @Test
    public void maximumFlow_gt() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector cmp = s.bv(4);
        Lit f = g.compareMaximumFlow(0, 3, cmp, Comparison.GT);
        assertEquals(s.solve(f), true);
        assertEquals(s.solve(f, cmp.eq(3)), false);
        assertEquals(s.solve(f, cmp.eq(2)), false);
        assertEquals(s.solve(f, cmp.eq(1)), true);

        assertEquals(s.solve(f, e_0_1.negate(), e_2_3.negate()), false);
        assertEquals(s.solve(f, e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(f, e_0_2.negate(), e_1_3.negate(), e_2_3.negate()), false);
        assertEquals(s.solve(f), true);
    }

    @Test
    public void maximumFlow_leq() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector cmp = s.bv(4);
        Lit f = g.compareMaximumFlow(0, 3, cmp, Comparison.LEQ);
        assertEquals(s.solve(f), true);
        assertEquals(s.solve(f, cmp.eq(3)), true);
        assertEquals(s.solve(f, cmp.eq(2)), true);
        assertEquals(s.solve(f, cmp.eq(1)), true);

        assertEquals(s.solve(f, e_0_1, e_0_2, e_1_3, e_2_3, cmp.eq(1)), false);
        assertEquals(s.solve(f, e_0_1, e_0_2, e_1_3, e_2_3, cmp.eq(2)), true);
        assertEquals(s.solve(f, e_0_1, e_0_2, e_1_3, e_2_3, cmp.eq(3)), true);
        assertEquals(s.solve(f, e_0_1.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(f, e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(f, e_0_2.negate(), e_1_3.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(f), true);
    }

    @Test
    public void maximumFlow_lt() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector cmp = s.bv(4);
        Lit f = g.compareMaximumFlow(0, 3, cmp, Comparison.LT);
        assertEquals(s.solve(f), true);
        assertEquals(s.solve(f, cmp.eq(3)), true);
        assertEquals(s.solve(f, cmp.eq(2)), true);
        assertEquals(s.solve(f, cmp.eq(1)), true);

        assertEquals(s.solve(f, e_0_1, e_0_2, e_1_3, e_2_3, cmp.eq(1)), false);
        assertEquals(s.solve(f, e_0_1, e_0_2, e_1_3, e_2_3, cmp.eq(2)), false);
        assertEquals(s.solve(f, e_0_1, e_0_2, e_1_3, e_2_3, cmp.eq(3)), true);

        assertEquals(s.solve(f, e_0_1.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(f, e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(f, e_0_2.negate(), e_1_3.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(f), true);
    }

    @Test
    public void maximumFlow() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));

        BitVector flow = g.maximumFlow(0, 3);
        assertEquals(s.solve(flow.gt(0)), true);
        assertEquals(s.solve(flow.eq(0)), true);
        assertEquals(s.solve(flow.eq(1)), true);
        assertEquals(s.solve(flow.eq(2)), true);
        assertEquals(s.solve(flow.eq(3)), false);

        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, flow.eq(1)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, flow.eq(2)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, flow.eq(3)), false);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate()), true);


    }


    @Test
    public void distance() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));

        BitVector dist = g.distance(0, 3);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), false);
        assertEquals(s.solve(dist.eq(1)), false);
        assertEquals(s.solve(dist.eq(2)), true);
        assertEquals(s.solve(dist.eq(3)), false);

        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), false);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), false);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate()), false);


    }

    @Test
    public void distanceConsts() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, 3);

        BitVector dist = g.distance(0, 3);
        assertEquals(s.solve(), true);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), false);
        assertEquals(s.solve(dist.eq(1)), false);
        assertEquals(s.solve(dist.eq(2)), true);
        assertEquals(s.solve(dist.eq(3)), false);

        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), false);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), false);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate()), false);


    }

    @Test
    public void distanceConstsGEQ() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, 3);

        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.GEQ);
        s.assertTrue(d);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), true);
        assertEquals(s.solve(dist.eq(1)), true);
        assertEquals(s.solve(dist.eq(2)), true);
        assertEquals(s.solve(dist.eq(3)), true);


        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), false);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)), true);

    }

    @Test
    public void distanceConstsGT() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, 3);

        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.GT);
        s.assertTrue(d);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), true);
        assertEquals(s.solve(dist.eq(1)), true);
        assertEquals(s.solve(dist.eq(2)), true);
        assertEquals(s.solve(dist.eq(3)), true);


        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), false);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)), true);

    }

    @Test
    public void distanceConstsLT() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, 3);

        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.LT);
        s.assertTrue(d);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), false);
        assertEquals(s.solve(dist.eq(1)), false);
        assertEquals(s.solve(dist.eq(2)), false);
        assertEquals(s.solve(dist.eq(3)), true);


        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), true);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), false);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)), false);

    }


    @Test
    public void distanceConstsLEQ() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, 3);

        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.LEQ);
        s.assertTrue(d);
        assertTrue(s.solve(dist.gt(0)));
        assertFalse(s.solve(dist.eq(0)));
        assertFalse(s.solve(dist.eq(1)));
        assertTrue(s.solve(dist.eq(2)));
        assertTrue(s.solve(dist.eq(3)));


        assertFalse(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)));
        assertTrue(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)));
        assertTrue(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)));

        assertFalse(s.solve(e_0_1.negate(), e_2_3.negate()));
        assertTrue(s.solve(e_0_2.negate(), e_1_3.negate()));
        assertFalse(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)));

    }

    @Test
    public void distanceGEQ() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.GEQ);
        s.assertTrue(d);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), true);
        assertEquals(s.solve(dist.eq(1)), true);
        assertEquals(s.solve(dist.eq(2)), true);
        assertEquals(s.solve(dist.eq(3)), true);


        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), false);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)), true);


    }

    @Test
    public void distanceGT() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.GT);
        s.assertTrue(d);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), true);
        assertEquals(s.solve(dist.eq(1)), true);
        assertEquals(s.solve(dist.eq(2)), true);
        assertEquals(s.solve(dist.eq(3)), true);


        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), false);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)), true);


    }

    @Test
    public void distanceLEQ() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.LEQ);
        s.assertTrue(d);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), false);
        assertEquals(s.solve(dist.eq(1)), false);
        assertEquals(s.solve(dist.eq(2)), true);
        assertEquals(s.solve(dist.eq(3)), true);


        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), true);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), true);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), false);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)), false);


    }


    @Test
    public void distanceLT() {
        Solver s = new Solver();

        Graph g = new Graph(s, 4);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph withone diagonal edge
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2, 2);
        Lit e_2_3 = g.addEdge(2, 3, s.bv(4, 1));
        BitVector dist = s.bv(4);
        Lit d = g.compareDistance(0, 3, dist, Comparison.LT);
        s.assertTrue(d);
        assertEquals(s.solve(dist.gt(0)), true);
        assertEquals(s.solve(dist.eq(0)), false);
        assertEquals(s.solve(dist.eq(1)), false);
        assertEquals(s.solve(dist.eq(2)), false);
        assertEquals(s.solve(dist.eq(3)), true);


        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)), false);
        assertEquals(s.solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)), true);

        assertEquals(s.solve(e_0_1.negate(), e_2_3.negate()), false);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate()), true);
        assertEquals(s.solve(e_0_2.negate(), e_1_3.negate(), e_2_3.negate(), dist.eq(1)), false);


    }


    @Test
    public void acyclicDirected() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }
        //create a directed square graph with two diagonal edges
        /*
         */
        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2);
        Lit e_2_3 = g.addEdge(2, 3);
        Lit e_3_0 = g.addEdge(3, 0);

        Lit r = g.acyclic();
        assertTrue(s.solve(r));
        assertTrue(s.solve(r, e_0_1.negate(), e_2_3.negate()));
        assertTrue(s.solve(r, e_0_2.negate(), e_1_3.negate()));
        assertTrue(s.solve(r, e_0_2.negate(), e_1_3.negate(), e_2_3.negate()));
        assertFalse(s.solve(r, e_0_1, e_1_2, e_2_3, e_3_0));
        assertTrue(s.solve(r));

        assertTrue(s.solve(r.negate()));
        assertFalse(s.solve(r.negate(), e_0_1.negate(), e_2_3.negate()));
        assertTrue(s.solve(r.negate(), e_0_2.negate(), e_1_3.negate()));
        assertFalse(s.solve(r.negate(), e_0_2.negate(), e_1_3.negate(), e_2_3.negate()));
        assertTrue(s.solve(r.negate(), e_0_1, e_1_2, e_2_3, e_3_0));
        assertTrue(s.solve(r.negate()));
        assertTrue(s.solve(r));

    }


    @Test
    public void acyclicUndirected() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for (int i = 0; i < 4; i++) {
            g.addNode();
        }

        //create a directed square graph with two diagonal edges
        /*
         */
        Lit e_0_1 = g.addEdge(0, 1);
        Lit e_0_2 = g.addEdge(0, 2);
        Lit e_1_3 = g.addEdge(1, 3);
        Lit e_1_2 = g.addEdge(1, 2);
        Lit e_2_3 = g.addEdge(2, 3);
        Lit e_3_0 = g.addEdge(3, 0);

        Lit r = g.acyclic();
        assertTrue(s.solve(r));
        assertTrue(s.solve(r, e_0_1.negate(), e_2_3.negate()));
        assertTrue(s.solve(r, e_0_2.negate(), e_1_3.negate()));
        assertTrue(s.solve(r, e_0_2.negate(), e_1_3.negate(), e_2_3.negate()));
        assertFalse(s.solve(r, e_0_1, e_1_2, e_2_3, e_3_0));
        assertTrue(s.solve(r));

        assertTrue(s.solve(r.negate()));
        assertFalse(s.solve(r.negate(), e_0_1.negate(), e_2_3.negate()));
        assertTrue(s.solve(r.negate(), e_0_2.negate(), e_1_3.negate()));
        assertFalse(s.solve(r.negate(), e_0_2.negate(), e_1_3.negate(), e_2_3.negate()));
        assertTrue(s.solve(r.negate(), e_0_1, e_1_2, e_2_3, e_3_0));
        assertTrue(s.solve(r.negate()));
        assertTrue(s.solve(r));

    }


    @Test
    public void edgeTest() {

    }
}