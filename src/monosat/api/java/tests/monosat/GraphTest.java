package monosat;

import org.junit.Test;

import static org.junit.Assert.*;

public class GraphTest {

    @Test
    public void nNodes() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for(int i = 0;i<10;i++){
            g.addNode();
            assertEquals(g.nNodes(),i+1);
        }
    }

    @Test
    public void nEdges() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for(int i = 0;i<4;i++){
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
        /*

        0 *--* 1
          |/ |
        2 *--* 3

         */
        Lit e_0_1 = g.addEdge(0,1);
        Lit e_0_2 = g.addEdge(0,2);
        Lit e_1_3 = g.addEdge(1,3);
        Lit e_1_2 = g.addEdge(1,2);
        Lit e_2_3 = g.addEdge(2,3);

        assertEquals(g.nEdges(),5);
    }

    @Test
    public void addNode() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for(int i = 0;i<10;i++){
            g.addNode();
            assertEquals(g.nNodes(),i+1);
        }
    }

    @Test
    public void bitwidth() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        assertEquals(g.bitwidth(),-1);
        Graph g2 = new Graph(s,4);
        assertEquals(g2.bitwidth(),4);
    }

    @Test
    public void addEdge() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for(int i = 0;i<4;i++){
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
        /*

        0 *--* 1
          |/ |
        2 *--* 3

         */
        Lit e_0_1 = g.addEdge(0,1);
        Lit e_0_2 = g.addEdge(0,2);
        Lit e_1_3 = g.addEdge(1,3);
        Lit e_1_2 = g.addEdge(1,2);
        Lit e_2_3 = g.addEdge(2,3);

        assertEquals(s.solve(e_0_1,e_0_2.negate(),e_1_3),true);
    }

    @Test
    public void reaches() {
        Solver s = new Solver();
        Graph g = new Graph(s);
        for(int i = 0;i<4;i++){
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
        /*

        0 *--* 1
          |/ |
        2 *--* 3

         */
        Lit e_0_1 = g.addEdge(0,1);
        Lit e_0_2 = g.addEdge(0,2);
        Lit e_1_3 = g.addEdge(1,3);
        Lit e_1_2 = g.addEdge(1,2);
        Lit e_2_3 = g.addEdge(2,3);

        Lit r = g.reaches(0,3);
        assertEquals(s.solve(r),true);
        assertEquals(s.solve(r,e_0_1.negate(),e_2_3.negate()),false);
        assertEquals(s.solve(r,e_0_2.negate(),e_1_3.negate()),true);
        assertEquals(s.solve(r,e_0_2.negate(),e_1_3.negate(),e_2_3.negate()),false);
        assertEquals(s.solve(r),true);
    }

    @Test
    public void maximumFlow_geq() {
          Solver s = new Solver();
          Graph g = new Graph(s,4);
          for(int i = 0;i<4;i++){
              g.addNode();
          }
          //create a directed square graph with one diagonal edges
          /*
  
          0 *--* 1
            |/ |
          2 *--* 3
  
           */
          
          Lit e_0_1 = g.addEdge(0,1);
          Lit e_0_2 = g.addEdge(0,2);
          Lit e_1_3 = g.addEdge(1,3);
          Lit e_1_2 = g.addEdge(1,2,2);
          Lit e_2_3 = g.addEdge(2,3,s.bv(4,1));
        BitVector cmp = s.bv(4);
        Lit f = g.maximumFlow_geq(0,3,cmp);
        assertEquals(s.solve(f),true);
        assertEquals(s.solve(f,cmp.eq(3)),false);
        assertEquals(s.solve(f,cmp.eq(2)),true);
        assertEquals(s.solve(f,cmp.eq(1)),true);
          assertEquals(s.solve(f,e_0_1.negate(),e_2_3.negate()),true);
        assertEquals(s.solve(f,e_0_1.negate(),e_2_3.negate(),cmp.eq(1)),false);
          assertEquals(s.solve(f,e_0_2.negate(),e_1_3.negate()),true);
          assertEquals(s.solve(f),true);
    }

    @Test
    public void maximumFlow_gt() {
        Solver s = new Solver();

        Graph g = new Graph(s,4);
        for(int i = 0;i<4;i++){
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0,1);
        Lit e_0_2 = g.addEdge(0,2);
        Lit e_1_3 = g.addEdge(1,3);
        Lit e_1_2 = g.addEdge(1,2,2);
        Lit e_2_3 = g.addEdge(2,3,s.bv(4,1));
        BitVector cmp = s.bv(4);
        Lit f = g.maximumFlow_gt(0,3,cmp);
        assertEquals(s.solve(f),true);
        assertEquals(s.solve(f,cmp.eq(3)),false);
        assertEquals(s.solve(f,cmp.eq(2)),false);
        assertEquals(s.solve(f,cmp.eq(1)),true);

        assertEquals(s.solve(f,e_0_1.negate(),e_2_3.negate()),false);
        assertEquals(s.solve(f,e_0_2.negate(),e_1_3.negate()),true);
        assertEquals(s.solve(f,e_0_2.negate(),e_1_3.negate(),e_2_3.negate()),false);
        assertEquals(s.solve(f),true);
    }

    @Test
    public void maximumFlow_leq() {
        Solver s = new Solver();

        Graph g = new Graph(s,4);
        for(int i = 0;i<4;i++){
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0,1);
        Lit e_0_2 = g.addEdge(0,2);
        Lit e_1_3 = g.addEdge(1,3);
        Lit e_1_2 = g.addEdge(1,2,2);
        Lit e_2_3 = g.addEdge(2,3,s.bv(4,1));
        BitVector cmp = s.bv(4);
        Lit f = g.maximumFlow_leq(0,3,cmp);
        assertEquals(s.solve(f),true);
        assertEquals(s.solve(f,cmp.eq(3)),true);
        assertEquals(s.solve(f,cmp.eq(2)),true);
        assertEquals(s.solve(f,cmp.eq(1)),true);

        assertEquals(s.solve(f,e_0_1,e_0_2,e_1_3,e_2_3 ,cmp.eq(1)),false);
        assertEquals(s.solve(f,e_0_1,e_0_2,e_1_3,e_2_3 ,cmp.eq(2)),true);
        assertEquals(s.solve(f,e_0_1,e_0_2,e_1_3,e_2_3 ,cmp.eq(3)),true);
        assertEquals(s.solve(f,e_0_1.negate(),e_2_3.negate()),true);
        assertEquals(s.solve(f,e_0_2.negate(),e_1_3.negate()),true);
        assertEquals(s.solve(f,e_0_2.negate(),e_1_3.negate(),e_2_3.negate()),true);
        assertEquals(s.solve(f),true);
    }

    @Test
    public void maximumFlow_lt() {
        Solver s = new Solver();

        Graph g = new Graph(s,4);
        for(int i = 0;i<4;i++){
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0,1);
        Lit e_0_2 = g.addEdge(0,2);
        Lit e_1_3 = g.addEdge(1,3);
        Lit e_1_2 = g.addEdge(1,2,2);
        Lit e_2_3 = g.addEdge(2,3,s.bv(4,1));
        BitVector cmp = s.bv(4);
        Lit f = g.maximumFlow_lt(0,3,cmp);
        assertEquals(s.solve(f),true);
        assertEquals(s.solve(f,cmp.eq(3)),true);
        assertEquals(s.solve(f,cmp.eq(2)),true);
        assertEquals(s.solve(f,cmp.eq(1)),true);

        assertEquals(s.solve(f,e_0_1,e_0_2,e_1_3,e_2_3 ,cmp.eq(1)),false);
        assertEquals(s.solve(f,e_0_1,e_0_2,e_1_3,e_2_3 ,cmp.eq(2)),false);
        assertEquals(s.solve(f,e_0_1,e_0_2,e_1_3,e_2_3 ,cmp.eq(3)),true);

        assertEquals(s.solve(f,e_0_1.negate(),e_2_3.negate()),true);
        assertEquals(s.solve(f,e_0_2.negate(),e_1_3.negate()),true);
        assertEquals(s.solve(f,e_0_2.negate(),e_1_3.negate(),e_2_3.negate()),true);
        assertEquals(s.solve(f),true);
    }

    @Test
    public void maximumFlow() {
        Solver s = new Solver();
        s.setOutputFile("/tmp/test.gnf");
        Graph g = new Graph(s,4);
        for(int i = 0;i<4;i++){
            g.addNode();
        }
        //create a directed square graph with one diagonal edges
          /*

          0 *--* 1
            |/ |
          2 *--* 3

           */

        Lit e_0_1 = g.addEdge(0,1);
        Lit e_0_2 = g.addEdge(0,2);
        Lit e_1_3 = g.addEdge(1,3);
        Lit e_1_2 = g.addEdge(1,2,2);
        Lit e_2_3 = g.addEdge(2,3,s.bv(4,1));

        BitVector flow = g.maximumFlow(0,3);
        assertEquals(s.solve(flow.gt(0)),true);
        assertEquals(s.solve(flow.eq(0)),true);
        assertEquals(s.solve(flow.eq(1)),true);
        assertEquals(s.solve(flow.eq(2)),true);
        assertEquals(s.solve(flow.eq(3)),false);

        assertEquals(s.solve(e_0_1,e_0_2,e_1_3,e_2_3 ,flow.eq(1)),false);
        assertEquals(s.solve(e_0_1,e_0_2,e_1_3,e_2_3 ,flow.eq(2)),true);
        assertEquals(s.solve(e_0_1,e_0_2,e_1_3,e_2_3 ,flow.eq(3)),false);

        assertEquals(s.solve(e_0_1.negate(),e_2_3.negate()),true);
        assertEquals(s.solve(e_0_2.negate(),e_1_3.negate()),true);
        assertEquals(s.solve(e_0_2.negate(),e_1_3.negate(),e_2_3.negate()),true);


    }
}