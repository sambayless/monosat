import jdk.nashorn.internal.runtime.regexp.joni.exception.ValueException;

/**
 * Literals are integers in the rance 0..nVars()*2
 * Literals come in pairs; positive literals are even, negative literals are odd.
 * Finally, there are two reserved literals:
 * lit_Undef and lit_Error, with values -2 and -1, respectively.
 *
 * In an ideal world, this would be _just_ an integer, or a type-checked value class of size 32 bits.
 * However, to deal with
 */
final public class Lit {
    //The value of this literal
    protected int l=-2;
    private final Lit neg; //every literal also has a pointer to its negation.

    public final static Lit Undef= new Lit(-2);
    public final static Lit Error= new Lit(-1);

    /**
     * Typically, users will create literals using Solver.newLit(), rather than constructing them directly.
     */
    protected Lit(){
        this.l = -2;
        this.neg = Error;
    }
    protected Lit(int variable){
        assert(variable>=0);
        this.l =variable*2;
        this.neg = new Lit(l+1,this);
    }
    private Lit(int l, Lit neg){
        this.l =l;
        this.neg = neg;
    }


    /**
     * Check whether the literal has a negation sign.
     * @return True if the literal is negative, false if it is positive.
     */
    public boolean sign(){
        return (l&1)==1;
    }

    public Lit negate(){
        return neg;
    }
    public Lit abs(){
        if(sign()){
            return neg;
        }else{
            return this;
        }
    }
    /**
     * Convert the literal into dimacs format
     * @return
     */
    public int toDimacs(){
        return ((l/2)+1) * (this.sign()? -1:1);
    }

    /**
     * Convert the literal into a variable.
     * @return
     */
    public int toVar(){return l/2;}
    /**
     * Throws an exception if this is not a legal literal (eg, if the index is <0)
     */
    public void validate(){
        if (l<0){
            throw new ValueException("Invalid literal");
        }
    }
}
