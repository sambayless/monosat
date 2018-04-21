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

import monosat.Lit;
import monosat.MonosatJNI;
import monosat.Solver;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

public class BitVector {
    protected int id;
    private int width;
    private Solver solver;
    private ArrayList<Lit> bits = new ArrayList<Lit>();

    public BitVector(Solver solver, ArrayList<Lit> bits) {
        this.solver = solver;
        id = MonosatJNI.newBitvector(solver.solverPtr, solver.bvPtr, solver.getLitBuffer(bits), bits.size());
        width = bits.size();
        for (Lit l : bits) {
            bits.add(l);
        }
    }

    public BitVector(Solver solver, int width, long constant) {
        this.solver = solver;
        assert (constant >= 0);
        assert (constant < (1L << width));
        id = MonosatJNI.newBitvector_const(solver.solverPtr, solver.bvPtr, width, constant);
        this.width = width;
        for (int i = 0; i < width; i++) {
            if ((constant & (1 << i)) == 1) {
                bits.add(Lit.True);
            } else {
                bits.add(Lit.False);
            }
        }
    }

    /**
     * Creates a new BitVector if the specified width.
     * If introduceLiterals is true, then introduces width literals,
     * else creates a bitvector that has no literals associated with it.
     *
     * @param solver
     * @param width
     */
    public BitVector(Solver solver, int width, boolean introduceLiterals) {
        this.solver = solver;
        if (width < 0) {
            throw new IllegalArgumentException("BitVector width must be >=0");
        }
        this.width = width;
        if (!introduceLiterals) {
            id = MonosatJNI.newBitvector_anon(solver.solverPtr, solver.bvPtr, width);
        } else {
            for (int i = 0; i < width; i++) {
                bits.add(new Lit(solver));
            }
            id = MonosatJNI.newBitvector(solver.solverPtr, solver.bvPtr, solver.getVarBuffer(bits, 0), bits.size());
        }
    }

    public BitVector(Solver solver, int width) {
        this(solver, width, true);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        BitVector bitVector = (BitVector) o;
        return id == bitVector.id &&
                Objects.equals(solver, bitVector.solver);
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, solver);
    }

    public Solver getSolver(){
        return solver;
    }

    public List<Lit> getBits() {
        return Collections.unmodifiableList(bits);
    }

    public Lit get(int bit) {
        return bits.get(bit);
    }

    public int width() {
        return width;
    }

    public int size() {
        return width();
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit gt(BitVector compareTo) {
        int l = MonosatJNI.newBVComparison_bv_gt(solver.solverPtr, solver.bvPtr, this.id, compareTo.id);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit geq(BitVector compareTo) {
        int l = MonosatJNI.newBVComparison_bv_geq(solver.solverPtr, solver.bvPtr, this.id, compareTo.id);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit lt(BitVector compareTo) {
        int l = MonosatJNI.newBVComparison_bv_lt(solver.solverPtr, solver.bvPtr, this.id, compareTo.id);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit leq(BitVector compareTo) {
        int l = MonosatJNI.newBVComparison_bv_leq(solver.solverPtr, solver.bvPtr, this.id, compareTo.id);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit neq(BitVector compareTo) {
        int l = MonosatJNI.newBVComparison_bv_neq(solver.solverPtr, solver.bvPtr, this.id, compareTo.id);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit eq(BitVector compareTo) {
        int l = MonosatJNI.newBVComparison_bv_eq(solver.solverPtr, solver.bvPtr, this.id, compareTo.id);
        return solver.toLit(l);
    }

    //Constant comparisons

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit gt(long compareTo) {
        int l = MonosatJNI.newBVComparison_const_gt(solver.solverPtr, solver.bvPtr, this.id, compareTo);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit geq(long compareTo) {
        int l = MonosatJNI.newBVComparison_const_geq(solver.solverPtr, solver.bvPtr, this.id, compareTo);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit lt(long compareTo) {
        int l = MonosatJNI.newBVComparison_const_lt(solver.solverPtr, solver.bvPtr, this.id, compareTo);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit leq(long compareTo) {
        int l = MonosatJNI.newBVComparison_const_leq(solver.solverPtr, solver.bvPtr, this.id, compareTo);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit neq(long compareTo) {
        int l = MonosatJNI.newBVComparison_const_neq(solver.solverPtr, solver.bvPtr, this.id, compareTo);
        return solver.toLit(l);
    }

    /**
     * Returns a literal which evaluates to true if this comparison holds, false otherwise.
     *
     * @param compareTo
     * @return
     */
    public Lit eq(long compareTo) {
        int l = MonosatJNI.newBVComparison_const_eq(solver.solverPtr, solver.bvPtr, this.id, compareTo);
        return solver.toLit(l);
    }

    /**
     * Creates a new bitvector consisting of the bits [this[0],..,this[size-1],append[0],..,append[append.size()-1]]
     * Does not introduce any new literals.
     *
     * @param append BitVector to append.
     * @return
     */
    public BitVector append(BitVector append) {
        int w = width() + append.width();
        BitVector result = new BitVector(solver, w);
        MonosatJNI.bv_concat(solver.solverPtr, solver.bvPtr, this.id, append.id, result.id);
        return result;
    }

    /**
     * Create a new bitvector consisting of the bits [this[lower],..,this[upper-1]]
     *
     * @param lower
     * @param upper
     * @return
     */
    public BitVector slice(int lower, int upper) {
        int w = upper - lower;
        assert (w >= 0);
        BitVector result = new BitVector(solver, w);
        MonosatJNI.bv_slice(solver.solverPtr, solver.bvPtr, this.id, lower, upper - 1, result.id);
        return result;
    }

    /**
     * Returns a Bitvector that represents the non-wrapping two's complement addition
     * of this and other. To prevent wrapping, the solver will enforce that a+b<(1<<width()).
     *
     * @param other
     * @return
     */
    public BitVector add(BitVector other) {
        return solver.add(this,other);
    }


    /**
     * Returns a Bitvector that represents the non-wrapping two's complement subtraction
     * of this and other. To prevent wrapping, the solver will enforce that a-b>=0.
     *
     * @param other
     * @return
     */
    public BitVector subtract(BitVector other) {
        return solver.subtract(this,other);
    }

    /**
     * Return the value of this bitvector from the solver.
     * Sometimes, a range of values may be determined by the solver to be satisfying.
     * If getMaximumValue is true, then largest value in that range will be returned,
     * otherwise, the smallest value is returned (this is relevant if optimization queries are being performed).
     *
     * @param getMaximumValue
     * @return
     */
    public long value(boolean getMaximumValue) {
        return solver.getValue(this);
    }

    /**
     * Return the value of this bitvector from the solver.
     * Sometimes, a range of values may be determined by the solver to be satisfying.
     * In this case, the smallest value is returned (this is relevant if optimization queries are being performed).
     *
     * @return
     */
    public long value() {
        return this.value(false);
    }
}
