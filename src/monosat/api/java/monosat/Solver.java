import monosat.MonosatJNI;

import java.util.ArrayList;

public class Solver {

    protected long solverPtr=0; //handle to the underlying monsoat solver instance.
    protected long bvPtr=0;
    //Instantiate a new solver
    public Solver(){
        solverPtr = MonosatJNI.newSolver();
        assert(solverPtr!=0);
        initBV();
    }
    public Solver(String args){
        solverPtr = MonosatJNI.newSolver(args);
        assert(solverPtr!=0);
        initBV();
    }
    public Solver(ArrayList<String> args){
        String arg = "";
        for(String s:args){
            arg+= s + " ";
        }
        solverPtr = MonosatJNI.newSolver(arg);
        assert(solverPtr!=0);
        initBV();
    }

    /**
     * Instantiate a bitvector theory solver
     */
    private void initBV(){
        assert(bvPtr==0);
        bvPtr = MonosatJNI.initBVTheory(solverPtr);
        assert(bvPtr!=0);
    }

    public void setOutputFile(String file){
        MonosatJNI.setOutputFile(solverPtr,file);
    }
    public Lit newLit(){
        return newLit(true);
    }
    public Lit newLit(boolean decisionVar){
        int var = MonosatJNI.newVar(solverPtr);
        assert(var>=0);
        return new Lit(var);
    }

    void releaseLiteral(Lit l){
        l.validate();
        int x = l.l;
        assert(x>=0);
    }

    void setDecisionLiteral(Lit l){
        MonosatJNI.setDecisionVar(l.toVar());
    }


}
