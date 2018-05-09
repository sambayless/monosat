package examples;

import monosat.*;

import java.util.Optional;

import static monosat.Logic.*;

/**
 * This is a brief introduction to MonoSAT's Java API,
 * which you can use to conveniently create and solve formulas with MonoSAT.
 * <p>
 * Before going any further, see the installation instructions for the Java library in [README].
 * You will need to compile with -DJAVA=ON, you will need to include monosat.jar on the classpath, and you will need to
 * set -Djava.library.path so that Java can find libmonosat.so/dylib.
 * <p>
 * Also, notice the two import statements above. The recommended way to use MonoSAT is to
 * statically import the monosat.Logic.*, so that standard logic operations (and, or, etc) will
 * be available directly.
 */
public class Tutorial {

    public static void main(String[] args) {

        //Create a new SAT solver
        Solver s = new Solver();

        /*
        When you create a solver, you can optionally pass some options to the underlying solver.
        In many (but not all) cases, the following option, which enables theory-based decision heuristics,
        is particularly good for performance:
        */
        s = new Solver("-decide-theories"); //create a new solver, store it in s, and allow the previous solver
        // (which no longer has any variables referring to it) to be garbage collected.

        //Create some literals in the solver
        Lit a = new Lit(s);
        Lit b = new Lit(s);

        //c is true if a is true or b is false, and false otherwise
        Lit c = or(a, not(b));

        //Add a unit clause to the solver, enforcing that c must be true in any satisfying solution
        assertTrue(c);

        //Solve the instance in MonoSAT, returning either true if the instance is SAT, or false if it is UNSAT.
        if (s.solve()) {
            System.out.println("SAT");
            //After a satisfiable call to solve(), you can query the assignments given by the solver to
            //individual variables using v.value()
            System.out.println("a: " + a.value());
            System.out.println("b: " + b.value());
            System.out.println("c: " + c.value());
        } else {
            System.out.println("UNSAT");
        }

        /*
        After a solve call, you can continue making further assertions, creating new variables, and making incremental
        calls to the solver.
        */

        Lit d= new Lit(s);
        assertTrue(implies(d, or(a,b)));

        /*
        There are also assertion forms for the common logic constructions, which are slightly more efficient than creating a
        new literal and asserting it to true, because they avoid introducing a new literal into the solver. A logically
        equivalent (but slightly more efficient) way to accomplish the above would have been:
        */
        assertImplies(d, or(a,b));

        /*
        Note that d does not yet have an assignment in the solver, and so calls to
        d.value() will throw an exception until the next solve call.
         */
        try{
            boolean r = d.value();
            System.out.println("This line is unreachable, because d.value() throws an exception");
        }catch(NoModelException e){

            //If you want to get the value of a literal that may or may not yet be assigned, use Lit.possibleValue();
            Optional<Boolean> val = d.possibleValue();
            System.out.println("Literal 'd' is not yet assigned, and so has value " + val);
        }


        if (s.solve()) {
            System.out.println("SAT");
            System.out.println("a: " + a.value());
            System.out.println("b: " + b.value());
            System.out.println("c: " + c.value());
            System.out.println("d: " + c.value()); //now d is assigned
        } else {
            System.out.println("UNSAT");
        }


        //Applying logic operators like 'and' and 'or' typically introduces a new literal and some constraints in
        //the solver. However, Lit.not/Logic.not is well optimized operator, and never introduces a new literal.
        assertTrue(not(and(a, b)));
        //The above is completely equivalent to
        assertTrue(nand(a, b));
        //And also completely equivalent to:
        assertTrue(not(not(not(and(a, b)))));
        // superfluous negations are eliminated by the solver and impose no overhead.

        if (s.solve()) {
            System.out.println("SAT");
            System.out.println("a: " + a.value());
            System.out.println("b: " + b.value());
            System.out.println("c: " + c.value());
            System.out.println("d: " + c.value()); //now d is assigned
        } else {
            System.out.println("UNSAT");
        }

        /*
        There is no way to remove assertions from MonoSAT yet, however, you can use assumptions to
        temporarily assert that a variable must be true (or false):
        */
        if (s.solve(b)) {
            System.out.println("SAT");
            System.out.println("a: " + a.value());
            System.out.println("b: " + b.value());
            System.out.println("c: " + c.value());
            System.out.println("d: " + c.value());
        } else {
            System.out.println("Temporarily UNSAT, under the assumption that 'b' is true");
        }

        /*
        If in the previous call, MonoSAT was only UNSAT under an assumption,
        then the solver may still be satisfiable in subsequent calls:
         */
        if (s.solve(not(b))) {

            System.out.println("SAT");
            System.out.println("a: " + a.value());
            System.out.println("b: " + b.value());
            System.out.println("c: " + c.value());
            System.out.println("d: " + c.value());
        } else {
            System.out.println("UNSAT");
        }

        // you can also assume multiple literals at once
        // the solver will then enforce the conjunction of these literals in the solver.
        if (s.solve(a,not(b), c)) {
            System.out.println("SAT");
            System.out.println("a: " + a.value());
            System.out.println("b: " + b.value());
            System.out.println("c: " + c.value());
            System.out.println("d: " + c.value());
        } else {
            System.out.println("UNSAT");
        }

        //***Theory Support***

        /*
        Now, onto the interesting stuff.
        In addition to Boolean logic, MonoSAT supports an extensive theory of finite graphs, including
        support for many common graph predicates such as reachability, shortest paths, maximum flows, and acyclicity.
        MonoSAT also has support for BitVectors and Cardinality/Pseudo-Boolean constraints.
        */

        // Constructing a graph in MonoSAT is as easy as:
        Graph g = new Graph(s);

        //Create three nodes
        int n0 = g.addNode();
        int n1 = g.addNode();
        int n2 = g.addNode();
        //nodes are just integers, with numbering starting from 0

        /*
        Add three directed edges to the graph.
        You can also create undirected edges, using g.addUndirectedEdge().
        */
        Lit e0 = g.addEdge(n0,n1);
        Lit e1 = g.addEdge(n1,n2);
        Lit e2 = g.addEdge(n0,n2);

        /*
         e0, e1, and e2 are *symbolic edges*, meaning that the edge (n0,n1) is included in G if and only if the
         theory atom (eg, the Literal) e0 is assigned to True by MonoSAT.
         You can use e0,e1, and e2 just like variables in MonoSAT, and in that way control which edges are in the graph
         using arbitrary Boolean logic.
         */

        assertNand(e0,e1,e2); //this is logically equivalent to assertTrue(not(and(e0,e1,e2)));
        assertOr(e0,e2);

        //You can even mix these symbolic edge variables with other logic from MonoSAT
        assertImplies(c, e0);

        /*
        Once you have created a graph and some edges, you can assert graph properties about that graph.
        For example, you can assert that node n2 must be reachable from node n0, in graph g
        */
        assertTrue(g.reaches(n0,n2));

        if (s.solve()) {
            System.out.println("SAT");
            System.out.println("e0: " + e0.value());
            System.out.println("e1: " + e1.value());
            System.out.println("e2: " + e2.value());
        }else {
            System.out.println("UNSAT");
        }

        //Graph predicates are 'double sided', so you can also assert that they are false, eg, in order to
        //prevent one node from reaching another:
        assertFalse(g.reaches(n1,n0));

        //Just like with edges, you can also mix graph predicates in with arbitrary logic anywhere you could use a Lit
        Lit r0_1 = g.reaches(n0,n1);
        assertOr(not(b), not(r0_1));

        if (s.solve()) {
            System.out.println("SAT");
            System.out.println("e0: " + e0.value());
            System.out.println("e1: " + e1.value());
            System.out.println("e2: " + e2.value());

            System.out.println("Use graph.draw() to draw the solution in graphviz/dot format:\n");
            System.out.println(g.draw());
        }else {
            System.out.println("UNSAT");
        }

        /*
        Edges can also have weights, represented as either constant longs,
        or as fixed-width, unsigned bounded bitvectors.

        (By unsigned bounded bitvectors, we mean that every bitvector in MonoSAT is asserted to
        be in the range [0, Max], and can never overflow/underflow -- if the only solution would involve overflowing,
         then the instance will instead be UNSAT)
        */
        //create a bitvector of width 4
        BitVector bv0 = new BitVector(s,4);
        BitVector bv1 = new BitVector(s,4);
        BitVector bv2 = new BitVector(s,4);

        /*
        BitVectors support addition, subtraction, and comparisons, but do not yet directly support
        negative values (the bitvectors are unsigned).
        */
        assertTrue(bv0.add(bv1).leq(7));
        assertTrue(bv0.add(bv2).geq( bv1));
        assertTrue(bv0.geq(2));

        if (s.solve()) {
            System.out.println("SAT");
            System.out.println("bv0: " + bv0.value());
            System.out.println("bv1: " + bv1.value());
            System.out.println("bv2: " + bv2.value());
        }else {
            System.out.println("UNSAT");
        }

        /*
        When creating an edge, you can use bitvectors or integers or longs
        as edge weights (otherwise, by default, every edge has constant weight '1'):
        */

        //Create a new graph with 3 nodes, that uses bitwidth 4 edges.
        Graph g2 = new Graph(s,4);
        int n4 = g2.addNode();
        int n5 = g2.addNode();
        int n6 = g2.addNode();

        //Add three weighted edges to the graph
        //Weights may be bitvectors, or integer constants.
        Lit e3 = g2.addEdge(n4,n5, bv0);
        Lit e4 = g2.addEdge(n5,n6, bv1);
        Lit e5 = g2.addEdge(n4,n6, bv2);

        /*
        Note that if you use bitvectors as edge weights, then all edges must use the same bitwidth, and you need to
        specify that bitwidth in the graph constructor. In this case, bv0, bv1, and bv2 are all bitwidth 4.
        */

        //MonoSAT supports several useful graph predicates in addition to reachability, including:

        /*
        Shortest path constraints:
        Assert that the distance from n0 to n2 is less or equal to 3
        */
        assertTrue(g2.compareDistance(n4,n6, LEQ,3));

        //You can also use BitVectors in the arguments of graph comparisons:
        BitVector bv3 = new BitVector(s,4);
        assertFalse(g2.compareDistance(n4,n6, LT,bv3));
        assertTrue(bv3.eq (bv0.add(bv1)));
        assertTrue(bv3.eq (2));

        /*
        You can also directly compute the distance as a bitvector which can then be used as normal
        (but be warned that this is substantially more expensive than a one-sided comparison)
        */
        BitVector dist = g2.distance(n4,n6);
        assertTrue(dist.eq(2));

        if (s.solve()) {
            System.out.println("SAT");
            System.out.println("e3: " + e3.value());
            System.out.println("e4: " + e4.value());
            System.out.println("e5: " + e5.value());
            System.out.println("bv0: " + bv0.value());
            System.out.println("bv1: " + bv1.value());
            System.out.println("bv2: " + bv2.value());
            System.out.println("bv3: " + bv3.value());
            System.out.println("distance: " + dist.value());
        }else {
            System.out.println("UNSAT");
        }

        /*
        MonoSAT also features highly optimized support for maximum flow constraints, allowing for comparisons against either a integer, or a bitvector:
        */
        assertTrue(g2.compareMaximumFlow(n4,n6,GEQ,3));

        BitVector bv4 = new BitVector(s,4);
        assertTrue(g2.compareMaximumFlow(n4,n6,GEQ,bv4));

        /*
        Just like with reachability and distance constraints, these maximum flow predicates are two sided
        so you can assert that the maximum flow must be less than a given bitvector, or you can include the
        maximum flow predicate as part of arbitrary Boolean logic.
        */
        assertOr(not(c),not(g2.compareMaximumFlow(n4,n6,GT,bv4.add(1))));

        /*
        And as with distance constraints,  you can also directly expose the maximum flow as a bitvector,
        (but this is substantially more expensive than simply performing a one sided comparison).
         */
        BitVector flow = g2.maximumFlow(n5,n6);

        if (s.solve()) {
            System.out.println("SAT");
            System.out.println("e3: " + e3.value());
            System.out.println("e4: " + e4.value());
            System.out.println("e5: " + e5.value());
            System.out.println("bv0: " + bv0.value());
            System.out.println("bv1: " + bv1.value());
            System.out.println("bv2: " + bv2.value());
            System.out.println("bv3: " + bv3.value());
            System.out.println("bv4: " + bv4.value());
            System.out.println("maximum flow: " + flow.value());
        }else {
            System.out.println("UNSAT");
        }



    }

}

