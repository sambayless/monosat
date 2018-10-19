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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

/**
 * A fixed-width BitVector representing an unsigned, non-wrapping integer value
 * in the range <code>0<=value<2^width<code/>.
 * <br>
 * Example usage:
 *
 * <blockquote><pre>{@code
 * Solver solver = new Solver();
 * BitVector a = new BitVector(solver,4); //create a new bitvector consisting of 4 literals.
 * BitVector b = a.add(1); //creates a new bitvector, equal to a+1.
 * BitVector c = new BitVector(solver,4);
 * BitVector d = a.subtract(c); //creates a new bitvector, equal to a-c.
 *
 * //Compare a to b and d
 * Lit x = a.lt(b);
 * Lit y = a.geq(d);
 * //Lit x will be true iff a < b, in the assignment chosen by the solver.
 * //Lit y will be true iff a >= d, in the assignment chosen by the solver.
 * //In this case, both comparisons are trivially true.
 *
 * Lit b0 = a.get(0); //Get the least significant bit of a.
 * Lit b3 = a.get(3); //Get the most significant bit of a.
 * }</pre></blockquote>
 *
 * BitVectors in MonoSAT are unsigned and non-wrapping, and logically equivalent
 * to mathematical integers that have been constrained to be in the range
 * 0<=value<2^width:
 * <br>
 * Asserting that a MonoSAT bitvector is outside of its legal range will result in an
 * unsatisfying formula.
 */
public final class BitVector {
  protected final int id;
  private final int width;
  private final Solver solver;
  private final ArrayList<Lit> bits = new ArrayList<Lit>();
  protected String _name;

  /**
   * MonoSAT supports different types of bitvectors.
   */
  public enum Type{
    /**
     * A theory-interpreted bitvector, with 1 literal bit bit
     */
    BV,
    /**
     * A bitvector that only operates at the SAT level until
     * it is completely assigned.
     */
     Lazy,
    /**
     * A bitvector that is theory-interpreted, but has no literals.
     */
    Anonymous
  }

  /**
   * Creates a new BitVector using the given literals. The width of the bitvector will equal
   * literals.size().
   *
   * @param solver The solver that this bitvector will belong to.
   * @param literals A non-empty set of at most 63 literals that will back this bitvector, in LSB
   *     order: literals.get(0) will represent the 1-bit of the bitvector, literals.get(1) the
   *     2-bit, etc.
   */
  public BitVector(Solver solver, ArrayList<Lit> literals) {
    this(solver, literals, "");
  }

  /**
   * Creates a new BitVector using the given literals. The width of the bitvector will equal
   * literals.size().
   *
   * @param solver The solver that this bitvector will belong to.
   * @param literals A non-empty set of at most 63 literals that will back this bitvector, in LSB
   *     order: literals.get(0) will represent the 1-bit of the bitvector, literals.get(1) the
   *     2-bit, etc.
   * @param name An (optional) name for the bitvector. Must be unique, and may not contain spaces or
   *     non-printable characters. May be the empty string (in which case the bitvector is unnamed).
   */
  public BitVector(Solver solver, ArrayList<Lit> literals, String name) {
    this.solver = solver;
    width = literals.size();
    if (width <= 0) {
      throw new IllegalArgumentException("BitVector must have a bit-width >= 0");
    } else if (width > 63) {
      throw new IllegalArgumentException("BitVector must have a bit-width <= 63");
    }
    id =
        MonosatJNI.newBitvector(
            solver.getSolverPtr(), solver.bvPtr, solver.getVarBuffer(literals,0), literals.size());
    //bits.addAll(literals);
    for (int i = 0; i < MonosatJNI.nBitvectorBits(solver.getSolverPtr(), solver.bvPtr, id); i++) {
      int l = MonosatJNI.getBitvectorBit(solver.getSolverPtr(), solver.bvPtr, id, i);
      bits.add(solver.getLiteral(l));
    }
    solver.registerBitVector(this);
    this._name = name;
    if (name.length() > 0) {
      MonosatJNI.setBitvectorName(solver.getSolverPtr(), solver.bvPtr, id, MonosatJNI.validID(name));
    }
  }

  /**
   * Creates a new BitVector of the specified width, with a constant value. If introduceLiterals is
   * true, then introduces width literals, else creates a bitvector that has no literals associated
   * with it.
   *
   * @param solver The solver that this bitvector will belong to.
   * @param width The number of bits in this BitVector. Width must be a non-zero positive integer <=
   *     63.
   * @param constant A non-negative constant value that this BitVector will represent. constant must
   *     be >=0, and < 1<<width.
   */
  public BitVector(Solver solver, int width, long constant) {
    this.solver = solver;
    if (width <= 0) {
      throw new IllegalArgumentException("BitVector must have a bit-width > 0");
    } else if (width > 63) {
      throw new IllegalArgumentException("BitVector must have a bit-width <= 63");
    }
    if (constant < 0) {
      throw new IllegalArgumentException("BitVectors can only represent values >=0");
    }
    if (constant >= (1L << width)) {
      throw new IllegalArgumentException("BitVectors can only represent values <= (2^width-1)");
    }
    id = MonosatJNI.newBitvector_const(solver.getSolverPtr(), solver.bvPtr, width, constant);
    this.width = width;
    for (int i = 0; i < MonosatJNI.nBitvectorBits(solver.getSolverPtr(), solver.bvPtr, id); i++) {
      int l = MonosatJNI.getBitvectorBit(solver.getSolverPtr(), solver.bvPtr, id, i);
      bits.add(solver.getLiteral(l));
    }
    /*   for (int i = 0; i < width; i++) {
      if ((constant & (1 << i)) == 1) {
        bits.add(Lit.True);
      } else {
        bits.add(Lit.False);
      }
    }*/
    this._name = "";
    solver.registerBitVector(this);
  }

  /**
   * Creates a new BitVector of the specified width. If introduceLiterals is true, then introduces
   * width literals, else creates a bitvector that has no literals associated with it.
   *
   * @param solver The solver that this bitvector will belong to.
   * @param width The number of bits in this BitVector. Width must be a non-zero positive integer <=
   *     63.
   * @param type The type of Bitvector to create.
   */
  public BitVector(Solver solver, int width, Type type) {
    this(solver, width, "", type);
  }

  /**
   * Creates a new BitVector of the specified width. If introduceLiterals is true, then introduces
   * width literals, else creates a bitvector that has no literals associated with it.
   *
   * @param solver The solver that this bitvector will belong to.
   * @param width The number of bits in this BitVector. Width must be a non-zero positive integer <=
   *     63.
   * @param name An (optional) name for the bitvector. Must be unique, and may not contain spaces or
   *     non-printable * characters. May be the empty string (in which case the bitvector is
   *     unnamed).
   */
  public BitVector(Solver solver, int width, String name) {
    this(solver, width, name, Type.BV);
  }

  /**
   * Creates a new BitVector of the specified width. If introduceLiterals is true, then introduces
   * width literals, else creates a bitvector that has no literals associated with it.
   *
   * @param solver The solver that this bitvector will belong to.
   * @param width The number of bits in this BitVector. Width must be a non-zero positive integer <=
   *     63.
   */
  public BitVector(Solver solver, int width) {
    this(solver, width, "", Type.BV);
  }

  /**
   * Creates a new BitVector of the specified width. If introduceLiterals is true, then introduces
   * width literals, else creates a bitvector that has no literals associated with it.
   *
   * @param solver The solver that this bitvector will belong to.
   * @param width The number of bits in this BitVector. Width must be a non-zero positive integer <=
   *     63.
   * @param name An (optional) name for the bitvector. Must be unique, and may not contain spaces or
   *     non-printable * characters. May be the empty string (in which case the bitvector is
   *     unnamed).
   * @param type The type of bitvector to create.
   */
  public BitVector(Solver solver, int width, String name, Type type) {
    this.solver = solver;
    if (width <= 0) {
      throw new IllegalArgumentException("BitVector width must be >0");
    } else if (width > 63) {
      throw new IllegalArgumentException("BitVector must have a bit-width <= 63");
    }
    this.width = width;
    if (name.length() > 0) {
      this._name = MonosatJNI.validID(name);

      if (MonosatJNI.hasBitvectorWithName(solver.getSolverPtr(), solver.bvPtr, name)) {
        // this name is already used
        throw new IllegalArgumentException("No two bitvectors may have the same (non-empty) name: " + name);
      }

    } else {
      this._name = "";
    }
    if (type==Type.Anonymous) {
      id = MonosatJNI.newBitvector_anon(solver.getSolverPtr(), solver.bvPtr, width);
    } else {
      for (int i = 0; i < width; i++) {
        bits.add(new Lit(solver));
      }
      if(type==Type.Lazy){
        id = MonosatJNI.newBitvector_lazy(solver.getSolverPtr(), solver.bvPtr, solver.getVarBuffer(bits, 0), bits.size
                ());
      }else if (type==Type.BV) {
        id = MonosatJNI.newBitvector(solver.getSolverPtr(), solver.bvPtr, solver.getVarBuffer(bits, 0), bits.size());
      }else{
        throw new IllegalArgumentException("No such type: " + type);
      }
    }
    assert (id >= 0);

    if (this._name.length() > 0) {
      MonosatJNI.setBitvectorName(solver.getSolverPtr(), solver.bvPtr, id, this._name);
    }
    solver.registerBitVector(this);
  }

  /**
   * Private constructor; create a Bitvector object to represent an already existing bitvector in
   * the solver. Creates a new BitVector of the specified width. If introduceLiterals is true, then
   * introduces width literals, else creates a bitvector that has no literals associated with it.
   *
   * @param solver The solver that this bitvector will belong to.
   * @param ignore This parameter exists only to distinguish this constructor from the publicly
   *     available ones. Ignore it.
   * @param bvID The ID of the bitvector in the solver
   */
  protected BitVector(Solver solver, Solver ignore, int bvID) {
    this.solver = solver;
    this.id = bvID;
    this.width = MonosatJNI.getBitvectorWidth(solver.getSolverPtr(), solver.bvPtr, bvID);
    for (int i = 0; i < MonosatJNI.nBitvectorBits(solver.getSolverPtr(), solver.bvPtr, bvID); i++) {
      int l = MonosatJNI.getBitvectorBit(solver.getSolverPtr(), solver.bvPtr, bvID, i);
      bits.add(solver.getLiteral(l));
    }
    _name = MonosatJNI.getBitvectorName(solver.getSolverPtr(), solver.bvPtr, id,0);
    solver.registerBitVector(this);
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    BitVector bitVector = (BitVector) o;
    return id == bitVector.id && Objects.equals(solver, bitVector.solver);
  }

  @Override
  public int hashCode() {
    return Objects.hash(id, solver);
  }

  /**
   * Get the solver that this BitVector belongs to.
   *
   * @return The solver that this BitVector belongs to.
   */
  public Solver getSolver() {
    return solver;
  }

  /**
   * Return an unmodifiable view of the literal that make up this bitvector (if any). The returned
   * list will either have exactly width literals, or be empty.
   *
   * @return The literals (if any) that make up this BitVector.
   */
  public List<Lit> getBits() {
    return Collections.unmodifiableList(bits);
  }

  /**
   * Get the literal at index 'bit' from the list of bits backing this bitvector.
   *
   * @param bit The index of the literal to retrieve.
   * @return The literal representing bit 'bit'.
   */
  public Lit get(int bit) {
    return bits.get(bit);
  }

  /**
   * Get the bitwidth of this (eg, number of bits) of this bitvector.
   *
   * @return The bitwidth of this bitvector.
   */
  public int width() {
    return width;
  }

  /**
   * Get the number of defined, explicit bit literals in this Bitvector.
   * Normally, this is equal to the width of the bitvector,
   * but in some cases a bitvector may have no explicit literals (in which case,
   * the size is exactly 0, even if the width of the bitvector is non-zero).
   * Size is always either equal to width(), or exactly 0.
   * @return The number of explicitly defined bit literals in the soover.
   */
  public int size() {
    return bits.size();
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is greater than compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = new BitVector(solver,4);
   * Lit l = a.gt(b);
   * //Lit l will be true iff a > b, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The BitVector that this BitVector will be compared to.
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit gt(BitVector compareTo) {
    return compare(Comparison.GT, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is greater or equal to compareTo,
   * false otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = new BitVector(solver,4);
   * Lit l = a.geq(b);
   * //Lit l will be true iff a >= b, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The BitVector that this BitVector will be compared to.
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit geq(BitVector compareTo) {
    return compare(Comparison.GEQ, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is less than compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = new BitVector(solver,4);
   * Lit l = a.lt(b);
   * //Lit l will be true iff a < b, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The BitVector that this BitVector will be compared to.
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit lt(BitVector compareTo) {
    return compare(Comparison.LT, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is less or equal to compareTo,
   * false otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = new BitVector(solver,4);
   * Lit l = a.leq(b);
   * //Lit l will be true iff a <= b, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The BitVector that this BitVector will be compared to.
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit leq(BitVector compareTo) {
    return compare(Comparison.LEQ, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is not equal to compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = new BitVector(solver,4);
   * Lit l = a.neq(b);
   * //Lit l will be true iff a != b, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The BitVector that this BitVector will be compared to.
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit neq(BitVector compareTo) {
    return compare(Comparison.NEQ, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is equal to compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = new BitVector(solver,4);
   * Lit l = a.eq(b);
   * //Lit l will be true iff a == b, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The BitVector that this BitVector will be compared to.
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit eq(BitVector compareTo) {
    return compare(Comparison.EQ, compareTo);
  }

  // Constant comparisons

  /**
   * Returns a literal which evaluates to true if this bitvector is greater than compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * Lit l = a.gt(5);
   * //Lit l will be true iff a > 5, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The non-negative long that this bitvector will be compared to. compareTo must
   *     be < 1<<width().
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit gt(long compareTo) {
    return compare(Comparison.GT, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is greater or equal to compareTo,
   * false otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * Lit l = a.geq(5);
   * //Lit l will be true iff a >= 5, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The non-negative long that this bitvector will be compared to. compareTo must
   *     be < 1<<width().
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit geq(long compareTo) {
    return compare(Comparison.GEQ, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is less than compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * Lit l = a.lt(5);
   * //Lit l will be true iff a < 5, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The non-negative long that this bitvector will be compared to. compareTo must
   *     be < 1<<width().
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit lt(long compareTo) {
    return compare(Comparison.LT, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is less or equal to compareTo,
   * false otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * Lit l = a.leq(5);
   * //Lit l will be true iff a <= 5, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The non-negative long that this bitvector will be compared to. compareTo must
   *     be < 1<<width().
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit leq(long compareTo) {
    return compare(Comparison.LEQ, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is not equal to compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * Lit l = a.neq(5);
   * //Lit l will be true iff a != 5, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The non-negative long that this bitvector will be compared to. compareTo must
   *     be < 1<<width().
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit neq(long compareTo) {
    return compare(Comparison.NEQ, compareTo);
  }

  /**
   * Returns a literal which evaluates to true if this bitvector is equal to compareTo, false
   * otherwise.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * Lit l = a.eq(5);
   * //Lit l will be true iff a == 5, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * @param compareTo The non-negative long that this bitvector will be compared to. compareTo must
   *     be < 1<<width().
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit eq(long compareTo) {
    return compare(Comparison.EQ, compareTo);
  }

  /**
   * Compare this bitvector to the bitvector 'compareTo'. Returns a literal that will evaluate to
   * true iff the comparison holds.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = new BitVector(solver,4);
   * Lit l = a.compare(Comparison.LT,b);
   * //Lit l will be true iff a < b, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * See also the short forms: BitVector.gt(),geq(),lt(),leq(),eq(),neq().
   *
   * @param c type of comparison to perform (GEQ,EQ,LT, etc.)
   * @param compareTo The BitVector that this BitVector will be compared to.
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit compare(Comparison c, BitVector compareTo) {
    switch (c) {
      case GT:
        return solver.toLit(
            MonosatJNI.newBVComparison_bv_gt(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo.id));
      case GEQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_bv_geq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo.id));
      case LT:
        return solver.toLit(
            MonosatJNI.newBVComparison_bv_lt(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo.id));
      case LEQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_bv_leq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo.id));
      case EQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_bv_eq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo.id));
      case NEQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_bv_neq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo.id));
    }
    // Comparison c should never be null
    throw new NullPointerException();
  }

  /**
   * Compare this bitvector to the bitvector 'compareTo'. Returns a literal that will evaluate to
   * true iff the comparison holds.
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * Lit l = a.compare(Comparison.EQ,5);
   * //Lit l will be true iff a == 5, in the assignment chosen by the solver.
   * }</pre>
   *
   * </blockquote>
   *
   * See also the short forms: BitVector.gt(),geq(),lt(),leq(),eq(),neq().
   *
   * @param c type of comparison to perform (GEQ,EQ,LT, etc.)
   * @param compareTo The non-negative long that this bitvector will be compared to. compareTo must
   *     be < 1<<width().
   * @return A literal that will evaluate to true iff the comparison holds.
   */
  public Lit compare(Comparison c, long compareTo) {
    switch (c) {
      case GT:
        return solver.toLit(
            MonosatJNI.newBVComparison_const_gt(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo));
      case GEQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_const_geq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo));
      case LT:
        return solver.toLit(
            MonosatJNI.newBVComparison_const_lt(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo));
      case LEQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_const_leq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo));
      case EQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_const_eq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo));
      case NEQ:
        return solver.toLit(
            MonosatJNI.newBVComparison_const_neq(
                solver.getSolverPtr(), solver.bvPtr, this.id, compareTo));
    }
    // Comparison c should never be null
    throw new NullPointerException();
  }

  /**
   * Creates a new bitvector consisting of the bits
   * [this[0],..,this[size-1],append[0],..,append[append.size()-1]] Does not introduce any new
   * literals.
   *
   * @param append BitVector to concatenate to this one.
   * @return A new BitVector, consisting of the concatenation of this bitvector and 'append'
   */
  public BitVector concatenate(BitVector append) {
    int w = width() + append.width();
    BitVector result = new BitVector(solver, w);
    MonosatJNI.bv_concat(solver.getSolverPtr(), solver.bvPtr, this.id, append.id, result.id);
    return result;
  }

  /**
   * Create a new BitVector consisting of the bits [this[lower],..,this[upper-1]]
   *
   * <blockquote>
   *
   * <pre>{@code
   * BitVector a = new BitVector(solver,4);
   * BitVector b = a.slice(0,3); //new bitvector of size 3, consisting of bits a[0],a[1], a[2]
   * BitVector c = a.slice(1,4); //new bitvector of size 3, consisting of bits a[1],a[2], a[3]
   * BitVector d = a.slice(2,3); //new bitvector of size 1, consisting of bit a[2]
   * }</pre>
   *
   * </blockquote>
   *
   * @param lower The (inclusive) lowest bit that will be the least significant bit of the new
   *     bitvector
   * @param upper The (exclusive) highest bit, one after the one that will be the most significant
   *     bit of the new bitvector
   * @return A new BitVector consisting of a sub-range of this bitvector
   */
  public BitVector slice(int lower, int upper) {
    int w = upper - lower;
    assert (w >= 0);
    BitVector result = new BitVector(solver, w);
    MonosatJNI.bv_slice(solver.getSolverPtr(), solver.bvPtr, this.id, lower, upper - 1, result.id);
    return result;
  }

  /**
   * Returns a Bitvector that represents the non-wrapping two's complement addition of this and
   * other. To prevent wrapping, the solver will enforce that a+b<2^width.
   *
   * @param other The bitvector to add to this one.
   * @return A Bitvector that represents this + other.
   */
  public BitVector add(BitVector other) {
    return solver.add(this, other);
  }

  /**
   * Returns a Bitvector that represents the non-wrapping two's complement addition of this and
   * other. To prevent wrapping, the solver will enforce that a+b<2^width.
   *
   * @param other The constant to add to this one.
   * @return A Bitvector that represents this + other.
   */
  public BitVector add(long other) {
    return solver.add(this, solver.bv(width(), other));
  }

  /**
   * Returns a Bitvector that represents the non-wrapping two's complement subtraction of this and
   * other. To prevent wrapping, the solver will enforce that a-b>=0.
   *
   * @param other The bitvector to subtract from this one.
   * @return A Bitvector that represents this - other.
   */
  public BitVector subtract(BitVector other) {
    return solver.subtract(this, other);
  }

  /**
   * Returns a Bitvector that represents the non-wrapping two's complement subtraction of this and
   * other. To prevent wrapping, the solver will enforce that a-b>=0.
   *
   * @param other The constant to subtract from this one.
   * @return A Bitvector that represents this - other.
   */
  public BitVector subtract(long other) {
    return solver.subtract(this, solver.bv(width(), other));
  }

  /**
   * Return the value of this bitvector from the solver. Sometimes, a range of values may be
   * determined by the solver to be satisfying. If getMaximumValue is true, then largest value in
   * that range will be returned, otherwise, the smallest value is returned (this is relevant if
   * optimization queries are being performed).
   *
   * @param getMaximumValue If true, and if the solver happened to that a range of values satisfied
   *     the formula, then return the largest value in that range (else, return the smallest).
   * @return A satisfying value for this bitvector, in the range 0<=value<(1<<width)
   */
  public long value(boolean getMaximumValue) {
    if (!MonosatJNI.hasModel(solver.getSolverPtr())) {
      throw new NoModelException(
          "Solver has no model (this may indicate either that the solve() has not yet been called, or that the most recent call to solve() returned a value other than true, or that a constraint was added into the solver after the last call to solve()).");
    }
    return MonosatJNI.getModel_BV(solver.getSolverPtr(), solver.bvPtr, id, getMaximumValue);
  }

  /**
   * Return the value of this bitvector from the solver. Sometimes, a range of values may be
   * determined by the solver to be satisfying. In this case, the smallest value is returned (this
   * is relevant if optimization queries are being performed).
   *
   * @return A satisfying value for this bitvector, in the range 0<=value<(1<<width)
   */
  public long value() {
    return this.value(false);
  }

  public String name() {
    return this._name;
  }

  public String toString() {
    if (this._name.length() > 0) {
      return _name + "(BV" + id + ",width" + width() + ")";
    } else {
      return "BV" + id + ",width" + width();
    }
  }

  /**
   * Get the integer ID of this bitvector
   *
   * @return the integer ID of this bitvector
   */
  protected int getID() {
    return this.id;
  }
}
