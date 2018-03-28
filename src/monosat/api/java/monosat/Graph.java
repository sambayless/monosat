package monosat;
import monosat.MonosatJNI;

public class Graph {
    Solver solver;
    private long graphPtr;
    int bitwidth = -1;

    public Graph(Solver solver){
        this.solver = solver;
        graphPtr = MonosatJNI.newGraph(solver.solverPtr);
    }
    public Graph(Solver solver, int bitwidth){
        this.solver = solver;
        graphPtr = MonosatJNI.newGraph(solver.solverPtr);
        assert(bitwidth>=0);
        this.bitwidth=bitwidth;
    }
    int nNodes(){
        return MonosatJNI.nNodes(solver.solverPtr,graphPtr);
    }
    int nEdges(){
        return MonosatJNI.nEdges(solver.solverPtr,graphPtr);
    }
    int addNode(){
        return MonosatJNI.newNode(solver.solverPtr,graphPtr);
    }

    int bitwidth(){
        return bitwidth;
    }

    Lit addEdge(int from, int to){
        return addEdge(from,to,1);
    }
    Lit addEdge(int from, int to, long constant_weight){
        if(bitwidth>=0) {
            return addEdge(from, to, solver.bv(bitwidth, constant_weight));
        }else{
            return solver.toLit(MonosatJNI.newEdge(solver.solverPtr,graphPtr,from,to,constant_weight));
        }
    }
    Lit addEdge(int from, int to, BitVector weight){
        if(this.bitwidth<0){
            throw new RuntimeException("In order to use bitvector edge weights, the bitwidth must be passed to the graph constructor, eg:" +
                    "Graph g = new Graph(solver, 8); will accept edges with bitvectors of size 8. Otherwise, edge weights are assumed to be constant integers.");
        }
        return solver.toLit(MonosatJNI.newEdge_bv(solver.solverPtr,graphPtr,from,to,weight.id));
    }
    Lit reaches(int from, int to){
        return solver.toLit(MonosatJNI.reaches(solver.solverPtr,graphPtr,from,to));
    }

    Lit compareDistance(int from, int to, long compareTo, Compare comparison){
        switch(comparison){
            case GEQ:
                return solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr,graphPtr,from,to,compareTo)).negate();
            case GT:
                return solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr,graphPtr,from,to,compareTo)).negate();
            case LEQ:
                return solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr,graphPtr,from,to,compareTo));
            case LT:
                return solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr,graphPtr,from,to,compareTo));
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.and(l1, l2);
            }case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_const(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_const(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }

    Lit compareDistance(int from, int to, BitVector compareTo, Compare comparison){
        switch(comparison){
            case GEQ:
                return solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr,graphPtr,from,to,compareTo.id)).negate();
            case GT:
                return solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr,graphPtr,from,to,compareTo.id)).negate();
            case LEQ:
                return solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr,graphPtr,from,to,compareTo.id));
            case LT:
                return solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr,graphPtr,from,to,compareTo.id));
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.and(l1, l2);
            }case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }

    BitVector distance(int from, int to){
        return distance(from,to,this.bitwidth);
    }
    BitVector distance(int from, int to, int bitwidth){
        if (bitwidth<0){
            //compute a 'large enough' bitwidth to hold the maximum possible distance
            bitwidth=8;//temporary
        }

        BitVector result = new BitVector(solver,bitwidth);

        Lit l1 = solver.toLit(MonosatJNI.shortestPath_leq_bv(solver.solverPtr, graphPtr, from, to, result.id));
        Lit l2 = solver.toLit(MonosatJNI.shortestPath_lt_bv(solver.solverPtr, graphPtr, from, to, result.id));
        //result is geq to the shortest path, and is not greater than the shortest path, and so it is exactly equal to the shortestPath
        solver.assertTrue(l1);
        solver.assertFalse(l2);
        return result;
    }


    /**
     * Compare a constant or bitvector to the maximum flow of a graph.
     * Note that if your goal is to assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     * @param from
     * @param to
     * @return
     */
    Lit compareMaximumFlow(int from, int to,long compareTo, Compare comparison){
        switch(comparison){
            case GEQ:
                return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr,graphPtr,from,to,compareTo));
            case GT:
                return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr,graphPtr,from,to,compareTo));
            case LEQ:
                return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr,graphPtr,from,to,compareTo)).negate();
            case LT:
                return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr,graphPtr,from,to,compareTo)).negate();
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.and(l1, l2);
            }case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr, graphPtr, from, to, compareTo));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr, graphPtr, from, to, compareTo));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }

    /**
     * Compare a constant or bitvector to the maximum flow of a graph.
     * Note that if your goal is to assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     * @param from
     * @param to
     * @return
     */
    Lit compareMaximumFlow(int from, int to,BitVector compareTo, Compare comparison){
        switch(comparison){
            case GEQ:
                return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr,graphPtr,from,to,compareTo.id));
            case GT:
                return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr,graphPtr,from,to,compareTo.id));
            case LEQ:
                return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr,graphPtr,from,to,compareTo.id)).negate();
            case LT:
                return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr,graphPtr,from,to,compareTo.id)).negate();
            case EQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.and(l1, l2);
            }case NEQ: {
                Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr, graphPtr, from, to, compareTo.id));
                return solver.nand(l1, l2);
            }
        }
        throw new RuntimeException("Unknown comparison");
    }
    BitVector maximumFlow(int from, int to){
        return maximumFlow(from,to,this.bitwidth);
    }
    /**
     * Create a new bitvector, equal in value to the maximum flow between from and to.
     * Note that if your goal is to  assert or reason about a one-sided comparison between the maximum flow and a constant or bitvector,
     * the compareMaximumFlow form is more efficient then the direct maximumFlow() form.
     * @param from
     * @param to
     * @return
     */
    BitVector maximumFlow(int from, int to, int bitwidth){
        if (bitwidth<0){
            //compute a 'large enough' bitwidth to hold the maximum possible flow
            bitwidth=8;//temporary
        }
        BitVector result = new BitVector(solver,bitwidth);
        Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr,graphPtr,from,to,result.id));
        Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr,graphPtr,from,to,result.id));
        //result is geq to the max flow, and is not greater than the max flow (and so it is exactly equal to the maxflow)
        solver.assertTrue(l1);
        solver.assertFalse(l2);
        return result;
    }


}
