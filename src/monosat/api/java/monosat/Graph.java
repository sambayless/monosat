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

    Lit maximumFlow_geq(int from, int to, BitVector compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr,graphPtr,from,to,compareTo.id));
    }
    Lit maximumFlow_gt(int from, int to, BitVector compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr,graphPtr,from,to,compareTo.id));
    }
    Lit maximumFlow_geq(int from, int to, long compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr,graphPtr,from,to,compareTo));
    }
    Lit maximumFlow_gt(int from, int to, long compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr,graphPtr,from,to,compareTo));
    }
    Lit maximumFlow_leq(int from, int to, BitVector compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr,graphPtr,from,to,compareTo.id)).negate();
    }
    Lit maximumFlow_lt(int from, int to, BitVector compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr,graphPtr,from,to,compareTo.id)).negate();
    }
    Lit maximumFlow_leq(int from, int to, long compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_gt(solver.solverPtr,graphPtr,from,to,compareTo)).negate();
    }
    Lit maximumFlow_lt(int from, int to, long compareTo){
        return solver.toLit(MonosatJNI.maximumFlow_geq(solver.solverPtr,graphPtr,from,to,compareTo)).negate();
    }
    /**
     * Create a new bitvector, equal in value to the maximum flow between from and to.
     * Note that if your goal is to assert a comparison between the maximum flow and a constant or bitvector,
     * one of the assertion forms will be much more efficient in the solver.
     * @param from
     * @param to
     * @return
     */
    BitVector maximumFlow(int from, int to){
        BitVector result = new BitVector(solver,bitwidth);
        Lit l1 = solver.toLit(MonosatJNI.maximumFlow_geq_bv(solver.solverPtr,graphPtr,from,to,result.id));
        Lit l2 = solver.toLit(MonosatJNI.maximumFlow_gt_bv(solver.solverPtr,graphPtr,from,to,result.id));
        Logic.AssertEqual(l1,l2);
        return result;
    }


}
