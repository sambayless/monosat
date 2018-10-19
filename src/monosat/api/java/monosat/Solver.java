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

import java.io.Closeable;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.*;
import java.util.logging.Logger;

/**
 * Represents a MonoSAT solver instance. Multiple solvers may be instantiated, with separate sets of
 * literals and clauses.
 *
 * <p>Example usage:
 *
 * <blockquote>
 *
 * <pre>{@code
 * Solver solver = new Solver();
 * Lit a = new Lit(solver);//create a new, free literal
 * BitVector b = new BitVector(solver,4);//create a BitVector of width 4.
 *
 * solver.addClause(a,b.gt(2)); //add a clause requiring either a to be true, or b to be greater than 2. *
 * solver.solve();//attempt to solve this formula
 *
 * //retrieve the values of the literal and bitvector from the solver's satisfying assignment
 * bool a_assign = a.value();
 * long b_assign = b.value();
 * }</pre>
 *
 * </blockquote>
 */
public final class Solver implements Closeable {
  /**
   * Holds weak references to all currently existing solvers, so that global logic operations on
   * True/False can be applied.
   */
  protected static final WeakHashMap<Solver, Boolean> solvers = new WeakHashMap<Solver, Boolean>();

  /** Logger for warnings in this library. */
  protected static final Logger log = Logger.getLogger("monosat");

  /** Largest constant BitVector value to cache. This may change in the future. */
  private static final long MAX_CACHE_CONST = 255;
  // Map from native pointers to graph objects, to avoid re-instantiating the same graph multiple
  // times
  protected final HashMap<Long, Graph> allGraphs = new HashMap<Long, Graph>();
  /** Caches instantiated, small BitVectors. */
  private final ArrayList<ArrayList<BitVector>> cached_bvs = new ArrayList<ArrayList<BitVector>>();
  /**
   * Holds instances of all literals, so that we don't need to create multiple literal objects for
   * the same literal
   */
  private final ArrayList<Lit> allLits = new ArrayList<>();
  /**
   * Contains only the positive versions of instantiated literals, in the order they were created.
   */
  private final ArrayList<Lit> positiveLiterals = new ArrayList<>();
  /** Each unique bitvector, stored by bvID. */
  private final ArrayList<BitVector> bvmap = new ArrayList<>();
  /**
   * Contains all unique BitVectors in this solver. No two bitvectors in this list have the same
   * bvID.
   */
  private final ArrayList<BitVector> allBVs = new ArrayList<>();
  /**
   * Handle to the underlying monosat solver instance. This is really a pointer, masquerading as a
   * long.
   */
  protected long solverPtr = 0;
  /**
   * Handle to the underlying monosat BitVector theory instance. This is really a pointer,
   * masquerading as a long.
   */
  protected long bvPtr = 0;
  /** Used internally to manage byte buffers for calls to the C api */
  private int buffer_size0 = 1024;
  /** Used internally to manage byte buffers for calls to the C api */
  private int buffer_size1 = 1024;
  /** Used internally to manage byte buffers for calls to the C api */
  private int buffer_size2 = 1024;
  /** Used internally to manage byte buffers for calls to the C api */
  private IntBuffer ints0;
  /** Used internally to manage byte buffers for calls to the C api */
  private IntBuffer ints1;
  /** Used internally to manage byte buffers for calls to the C api */
  private IntBuffer ints2;

  /** Instantiate a new Solver, with default settings. */
  public Solver() {
    this(false);
  }

  /**
   * Instantiate a new Solver, with the given settings. MonoSAT has many command line options, which
   * can be passed here. A typical usage is:
   *
   * <blockquote>
   *
   * <pre>{@code
   * Solver s = new Solver("-decide-theories"); //Create a solver with fast, graph-based decision heuristics enabled.
   * }</pre>
   *
   * </blockquote>
   */
  public Solver(String args) {
    this(args, false);
  }

  /**
   * Instantiate a new Solver, with default settings.
   *
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   */
  public Solver(boolean enablePreprocessing) {
    this("", enablePreprocessing);
  }
  /**
   * Instantiate a new Solver, with the given list of settings. MonoSAT has many command line
   * options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   */
  public Solver(String[] args) {
    this(Arrays.asList(args), false);
  }

  /**
   * Instantiate a new Solver, with the given list of settings, and with constraints written to the
   * specified output file. MonoSAT has many command line options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   * @param outputFile A file to write constraints to, in MonoSAT's GNF format.
   */
  public Solver(String[] args, String outputFile) {
    this(Arrays.asList(args), false, outputFile);
  }

  /**
   * Instantiate a new Solver, with the given list of settings. MonoSAT has many command line
   * options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   */
  public Solver(String[] args, boolean enablePreprocessing) {
    this(Arrays.asList(args), enablePreprocessing, "");
  }

  /**
   * Instantiate a new Solver, with the given list of settings, and with constraints written to the
   * specified output file. MonoSAT has many command line options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   * @param outputFile A file to write constraints to, in MonoSAT's GNF format.
   */
  public Solver(String[] args, boolean enablePreprocessing, String outputFile) {
    this(Arrays.asList(args), enablePreprocessing, outputFile);
  }

  /**
   * Instantiate a new Solver, with the given list of settings. MonoSAT has many command line
   * options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   */
  public Solver(List<String> args) {
    this(String.join(" ",args), false);
  }

  /**
   * Instantiate a new Solver, with the given list of settings, and with constraints written to the
   * specified output file. MonoSAT has many command line options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   * @param outputFile A file to write constraints to, in MonoSAT's GNF format.
   */
  public Solver(List<String> args, String outputFile) {
    this(String.join(" ",args), false, outputFile);
  }

  /**
   * Instantiate a new Solver, with the given list of settings. MonoSAT has many command line
   * options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   */
  public Solver(List<String> args, boolean enablePreprocessing) {
    this(String.join(" ",args), enablePreprocessing, "");
  }

  /**
   * Instantiate a new Solver, with the given list of settings, and with constraints written to the
   * specified output file. MonoSAT has many command line options, which can be passed here.
   *
   * @param args A list of strings, with each string representing an argument to pass to the solver.
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   * @param outputFile A file to write constraints to, in MonoSAT's GNF format.
   */
  public Solver(List<String> args, boolean enablePreprocessing, String outputFile) {
    this(String.join(" ",args), enablePreprocessing, outputFile);
  }

  /**
   * Instantiate a new Solver, with the given settings. MonoSAT has many command line options, which
   * can be passed here.
   *
   * @param args A string of space-delimited arguments to pass to the solver.
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   */
  public Solver(String args, boolean enablePreprocessing) {
    this(args, enablePreprocessing, "");
  }

  /**
   * Instantiate a new Solver, with the given settings. MonoSAT has many command line options, which
   * can be passed here.
   *
   * @param args A string of space-delimited arguments to pass to the solver.
   * @param outputFile A file to write constraints to, in MonoSAT's GNF format.
   */
  public Solver(String args, String outputFile) {
    this(args, false, outputFile);
  }
  private String outputFile = "";
  /**
   * Instantiate a new Solver, with the given settings. MonoSAT has many command line options, which
   * can be passed here.
   *
   * @param args A string of space-delimited arguments to pass to the solver.
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   * @param outputFile A file to write constraints to, in MonoSAT's GNF format.
   */
  public Solver(String args, boolean enablePreprocessing, String outputFile) {
    this(args,enablePreprocessing, outputFile,true);
    }

  /**
   * Instantiate a new Solver, with the given settings. MonoSAT has many command line options, which
   * can be passed here.
   *
   * @param args A string of space-delimited arguments to pass to the solver.
   * @param enablePreprocessing If true, support for pre-processing is enabled. Pre-processing can
   *     improve performance, but requires some extra care in solver usage, and so is disabled by
   *     default.
   * @param outputFile A file to write constraints to, in MonoSAT's GNF format.
   * @param defineTrue If true (the default), define the global literal 'True' in the solver. If false,
   * then the global literal True will be left unconstrained in the solver. This is
   * occasionally required for dealing with plain CNFs, in which '1' is not defined to be True.
   *
   */
  public Solver(String args, boolean enablePreprocessing, String outputFile, boolean defineTrue) {
    if (args != null && args.length() > 0) {
      solverPtr = MonosatJNI.newSolver("monosat " + args);
    } else {
      solverPtr = MonosatJNI.newSolver();
    }
    if (solverPtr == 0) {
      throw new RuntimeException("Failed to created solver");
    }
    // Keep a global list of all solvers that yet created
    solvers.put(this, true);
    this.outputFile= outputFile;
    if (outputFile != null && outputFile.length() > 0) {
      MonosatJNI.setOutputFile(solverPtr, outputFile);
    }
    if (!enablePreprocessing) {
      disablePreprocessing();
    }
    if (defineTrue) {
      int true_lit = MonosatJNI.true_lit(solverPtr);
      assert (true_lit == 0);
      registerLit(Lit.True, Lit.False);
      this.addClause(Lit.True);
      // Ensure that the True lit returned by circuit operations (eg, and())
      // is the same as this True Lit.

      assert (Lit.True == toLit(true_lit));

    }
    initBV();
    initBuffers();


  }

    /**
     * Get the native pointer to this solver, or throw an NPE
     * @return The native pointer for this solver object, if it exists.
     * @throws NullPointerException If the solver has been deleted.
     */
  protected long getSolverPtr(){
      if(solverPtr==0){
          throw new NullPointerException("Native solver object has been deleted.");
      }else{
          return solverPtr;
      }
  }

  /** @return The version string of the MonoSAT native library. */
  public static String getVersion() {
    return MonosatJNI.getVersion();
  }



  /**
   * Load a formula into the solver. The formula may be in standard DIMACS CNF format, or MonoSAT's
   * extended GNF format.
   *
   * <p>The loaded constraints will exist alongside any constraints already in the solver (eg, any
   * previously defined literals will still be defined, and if the solver was previously UNSAT, it
   * will remain UNSAT after loading new constraints).
   *
   * @param filename The CNF or GNF constraints to read.
   * @throws IOException If the file could not be read.
   */
  public void loadConstraints(String filename) throws IOException {
      //note that any embedded solve calls in the GNF file will be ignored
    MonosatJNI.loadGNF(getSolverPtr(), filename);
  }

  /** If the solver is writing constraints to a file, close that file. */
  public void closeConstraintFile() {
    MonosatJNI.closeFile(getSolverPtr());
  }

  /** If the solver is writing constraints to a file, flush that file. */
  public void flushConstraintFile() {
    MonosatJNI.flushFile(getSolverPtr());
  }

  /**
   * Release all native resources associated with this solver. This does not normally need to be
   * called manually (though it is safe to do so), as it will be called automatically during garbage
   * collection.
   *
   * It is safe to close a solver multiple times.
   *
   * Note that many methods of the solver (and of objects belonging to the solver, such as literals or bitvectors or
   * graphs) may throw NullPointerExceptions if accessed after closing the solver.
   */
  @Override
  public synchronized void close() {
    // Does this method actually need to be synchronized?
    //It is always safe to call this method multiple times
    if (solverPtr != 0) {
      solvers.remove(this);
      MonosatJNI.deleteSolver(getSolverPtr());
      solverPtr = 0;
    }
  }

  @Override
  protected void finalize() {
    // Finalize is problematic, but we follow the Java standard library's
    // FileInputStream's here to release native resources on finalize().

    // Consider replacing with Java 9's java.lang.ref.Cleaner in the future.
    // But for now, sticking with finalize (to support Java 8)
    if (solverPtr != 0) {
      close();
    }
  }

  /** Internal method to initialize buffers for passing arguments to and from the native library. */
  private void initBuffers() {
    initBuffer(0);
    initBuffer(1);
    initBuffer(2);
  }

  /**
   * Internal method to initialize buffers for passing arguments to and from the native library.
   *
   * @param bufferN A buffer index to initialize. Must be 0,1, or 2.
   */
  private void initBuffer(int bufferN) {
    if (bufferN == 0) {
      ByteBuffer b = ByteBuffer.allocateDirect(buffer_size0 * 4); // 4 bytes per integer
      b.order(ByteOrder.LITTLE_ENDIAN);
      ints0 = b.asIntBuffer();
    } else if (bufferN == 1) {
      ByteBuffer b = ByteBuffer.allocateDirect(buffer_size0 * 4);
      b.order(ByteOrder.LITTLE_ENDIAN);
      ints1 = b.asIntBuffer();
    } else if (bufferN == 2) {
      ByteBuffer b = ByteBuffer.allocateDirect(buffer_size0 * 4);
      b.order(ByteOrder.LITTLE_ENDIAN);
      ints2 = b.asIntBuffer();
    } else {
      throw new IllegalArgumentException("Unknown buffer " + bufferN);
    }
  }

  /**
   * Get an integer direct buffer for passing arguments to and from the native library.
   *
   * @param bufferN The index of the buffer (must be 0,1, or 2).
   * @param minsize The minimum size that the buffer must have allocated (must be >=0).
   * @return A direct integer buffer of sufficient size to hold minsize integers.
   */
  protected IntBuffer getBuffer(int bufferN, int minsize) {
    if (minsize < 0) {
      throw new IllegalArgumentException("minsize must be >=0");
    }
    if (bufferN == 0) {
      if (minsize >= buffer_size0) {
        buffer_size0 = minsize * 2;
        initBuffer(bufferN);
      }
      return ints0;
    } else if (bufferN == 1) {
      if (minsize >= buffer_size1) {
        buffer_size1 = minsize * 2;
        initBuffer(bufferN);
      }
      return ints1;
    } else if (bufferN == 2) {
      if (minsize >= buffer_size2) {
        buffer_size2 = minsize * 2;
        initBuffer(bufferN);
      }
      return ints2;
    }
    throw new IllegalArgumentException("BufferN must be between 0 and 2");
  }

  /** Instantiate the bitvector theory in this solver. */
  private void initBV() {
    assert (bvPtr == 0);
    bvPtr = MonosatJNI.initBVTheory(getSolverPtr());
    assert (bvPtr != 0);
  }

  /**
   * Throw an exception if this literal is invalid. Otherwise, do nothing.
   *
   * @param l A literal to check.
   * @throws NullPointerException If l is null.
   * @throws IllegalArgumentException If l does not belong to this solver or is not a valid literal.
   */
  protected void validate(Lit l) {
    if (l == Lit.True || l == Lit.False) return;
    if (l == null) {
      throw new NullPointerException("Literal is null");
    } else if (l.l < 0) {
      throw new IllegalArgumentException("Literal " + l.toString() + " is not a valid literal.");
    } else if (l.solver != this) {
      throw new IllegalArgumentException(
          "Cannot pass literal belonging to solver "
              + (l.solver == null ? "null" : l.solver.toString())
              + " to solver "
              + toString());
    } else if (l.toVar() >= nVars()) {
      throw new IllegalArgumentException("Literal is undefined in solver (too large)");
    }
  }

  /**
   * Check an array of literal arguments for invalid literals.
   *
   * @param args An array of literals to check.
   * @throws NullPointerException If any literal in args is null.
   * @throws IllegalArgumentException If any literal in args does not belong to this solver or is
   *     not a valid literal.
   */
  protected void validate(Lit... args) {
    for (Lit l : args) {
      validate(l);
    }
  }



  /**
   * Check a collection of literal arguments for invalid literals.
   *
   * @param args Literals to check.
   * @throws NullPointerException If any literal in args is null.
   * @throws IllegalArgumentException If any literal in args does not belong to this solver or is
   *     not a valid literal.
   */
  protected void validate(Iterable<Lit> args) {
    for (Lit l : args) {
      validate(l);
    }
  }

  /**
   * Check an array of BitVector arguments.
   *
   * @param args An array of BitVectors to check.
   * @throws NullPointerException If any BitVector in args is null.
   * @throws IllegalArgumentException If any BitVector in args does not belong to this solver or is
   *     not a valid BitVector.
   */
  protected void validate(BitVector... args) {
    for (BitVector bv : args) {
      if (bv == null) {
        throw new NullPointerException("BitVector is null");
      } else if (bv.getSolver() != this) {
        throw new IllegalArgumentException(
            "Cannot pass BitVector belonging to solver "
                + (bv.getSolver() == null ? "null" : bv.getSolver().toString())
                + " to solver "
                + toString());
      } else if (bv.id < 0) {
        throw new IllegalArgumentException("BitVector is undefined " + bv.toString());
      }
    }
  }

  /**
   * Check a collection of BitVector arguments.
   *
   * @param args A collection of BitVectors to check.
   * @throws NullPointerException If any BitVector in args is null.
   * @throws IllegalArgumentException If any BitVector in args does not belong to this solver or is
   *     not a valid BitVector.
   */
  protected void validateBV(Collection<BitVector> args) {
    for (BitVector bv : args) {
      if (bv == null) {
        throw new IllegalArgumentException("Literal is null");
      } else if (bv.getSolver() != this) {
        throw new IllegalArgumentException(
            "Cannot pass literal belonging to solver "
                + (bv.getSolver() == null ? "null" : bv.getSolver().toString())
                + " to solver "
                + toString());
      } else if (bv.id < 0) {
        throw new IllegalArgumentException("BitVector is undefined " + bv.toString());
      }
    }
  }

  /**
   * Release this literal in the solver, so that it can be re-used. Use this method with care, as l
   * may no longer be valid after this call!
   *
   * @param l The literal to release.
   */
  public void releaseLiteral(Lit l) {
    validate(l);
    int x = l.l;
    assert (x >= 0);
    MonosatJNI.releaseLiteral(getSolverPtr(), x);
  }

  /**
   * Set whether the solver can make decisions on the value of l or not. By default, all literals
   * are decidable.
   *
   * @param l The literal to set.
   * @param decidable If true, the literal can be used in decisions in the solver. else, it cannot
   *     be.
   */
  public void setDecisionLiteral(Lit l, boolean decidable) {
    validate(l);
    MonosatJNI.setDecisionVar(getSolverPtr(), l.toVar(), decidable);
  }

  /**
   * Check if a literal is decidable.
   *
   * @param l The literal to check.
   * @return True if l is decidable.
   */
  public boolean isDecisionLiteral(Lit l) {
    validate(l);
    return MonosatJNI.isDecisionVar(getSolverPtr(), l.toVar());
  }

  /**
   * Set the decision priority of a literal in the solver's decision heuristic.
   *
   * @param l The literal to check.
   * @param priority The integer decision priority.
   */
  public void setDecisionPriority(Lit l, int priority) {
    validate(l);
    MonosatJNI.setDecisionPriority(getSolverPtr(), l.toVar(), priority);
  }

  /**
   * Get the decision priority of a literal.
   *
   * @param l The literal whose priority should be returned.
   * @return The decision priority of literal 'l'.
   */
  public int getDecisionPriority(Lit l) {
    validate(l);
    return MonosatJNI.getDecisionPriority(getSolverPtr(), l.toVar());
  }

  /**
   * Set the decision polarity of a literal in the solver's decision heuristic.
   *
   * @param l The literal to check.
   * @param polarity The integer decision polarity.
   */
  public void setDecisionPolarity(Lit l, boolean polarity) {
    validate(l);
    MonosatJNI.setDecisionPolarity(getSolverPtr(), l.toVar(), polarity);
  }

  /**
   * Get the decision polarity of a literal.
   *
   * @param l The literal whose polarity should be returned.
   * @return The decision polarity of literal 'l'.
   */
  public boolean getDecisionPolarity(Lit l) {
    validate(l);
    return MonosatJNI.getDecisionPolarity(getSolverPtr(), l.toVar());
  }

  /**
   * If preprocessing is enabled in the solver (by default, it is not), then prevent this literal
   * from being simplified by preprocessing.
   *
   * @param l The literal to disable simplification on.
   * @throws RuntimeException If the Literal has already been eliminated by preprocessing.
   */
  public void disallowSimplification(Lit l) {
    validate(l);
    if (!MonosatJNI.disallowLiteralSimplification(getSolverPtr(), l.toVar())) {
      throw new RuntimeException("Literal has already been eliminated.");
    }
  }

  /**
   * If preprocessing was enabled in this solver, disable it. Otherwise, do nothing. Using
   * preprocessing requires some care, and so by default preprocessing is disabled.
   */
  public void disablePreprocessing() {
    MonosatJNI.disablePreprocessing(getSolverPtr());
  }

  /**
   * Get the number of variables in this solver. Note that variable 'True' is always defined in a
   * solver, so nVars() is always at least 1, even if no Literals have been explicitly created.
   *
   * @return The number of variables in this solver.
   */
  public int nVars() {
    return MonosatJNI.nVars(getSolverPtr());
  }

  /**
   * Get the number of (non-learned) clauses in this solver.
   *
   * @return The number of clauses in this solver.
   */
  public int nClauses() {
    return MonosatJNI.nClauses(getSolverPtr());
  }

  /**
   * Get the number of learned clauses in this solver.
   *
   * @return The number of learned clauses in this solver.
   */
  public int nLearnedClauses() {
    return MonosatJNI.nLearnedClauses(getSolverPtr());
  }

  /**
   * Get the number of bitvectors in this solver.
   *
   * @return The number of bitvectors in this solver.
   */
  public int nBitVectors() {
    return MonosatJNI.nBitvectors(getSolverPtr(), bvPtr);
  }

  /**
   * Internal method for converting java arrays of bitvectors into byte buffers (as required by the
   * Monosat JNI)
   *
   * @param bvs A collection of BitVectors whose IDs will be stored in the specified buffer.
   * @param bufferN The buffer to store BV IDs in.
   * @return The buffer that the BV ids were stored in.
   */
  protected IntBuffer getBVBuffer(Collection<BitVector> bvs, int bufferN) {
    assert (bufferN < 3);
    assert (bufferN >= 0);
    IntBuffer buffer = getBuffer(bufferN, bvs.size());
    int index = 0;
    for (BitVector bv : bvs) {

      buffer.put(index, bv.id);
      index++;
    }
    return buffer;
  }

  /**
   * Internal method for converting java collections of literals into byte buffers (as required by
   * the Monosat JNI)
   *
   * @param clause A collection of literals to store in a buffer.
   * @param bufferN The buffer to store the literals in.
   * @return A buffer containing the Variable representation of the specified literals.
   */
  protected IntBuffer getVarBuffer(Collection<Lit> clause, int bufferN) {
    assert (bufferN < 3);
    assert (bufferN >= 0);
    IntBuffer buffer = getBuffer(bufferN, clause.size());
    int index = 0;
    for (Lit l : clause) {
      buffer.put(index, l.toVar());
      index++;
    }
    return buffer;
  }

  /**
   * Internal method for converting java collections of Lits into byte buffers (as required by the
   * Monosat JNI)
   *
   * @param clause A collection of literals to store in a buffer.
   * @return A buffer containing the integer representation of the specified literals.
   */
  protected IntBuffer getLitBuffer(Collection<Lit> clause) {
    return getLitBuffer(clause, 0);
  }

  /**
   * Internal method for converting an array of Lits into byte buffers (as required by the Monosat
   * JNI)
   *
   * @param clause An array of literals to store in a buffer.
   * @param bufferN The index of the buffer to fill.
   * @return A buffer containing the integer representation of the specified literals.
   */
  protected IntBuffer getLitBuffer(Lit[] clause, int bufferN) {
    assert (bufferN < 3);
    assert (bufferN >= 0);
    IntBuffer buffer = getBuffer(bufferN, clause.length);
    int index = 0;
    for (Lit l : clause) {
      buffer.put(index, l.l);
      index++;
    }
    return buffer;
  }

  /**
   * Internal method for converting java collections of Lits into byte buffers (as required by the
   * Monosat JNI)
   *
   * @param clause An collection of literals to store in a buffer.
   * @param bufferN The index of the buffer to fill.
   * @return A buffer containing the integer representation of the specified literals.
   */
  protected IntBuffer getLitBuffer(Collection<Lit> clause, int bufferN) {
    assert (bufferN < 3);
    assert (bufferN >= 0);
    IntBuffer buffer = getBuffer(bufferN, clause.size());
    int index = 0;
    for (Lit l : clause) {
      buffer.put(index, l.l);
      index++;
    }
    return buffer;
  }

  /**
   * Internal method for converting java collections of integers into byte buffers (as required by
   * the Monosat JNI)
   *
   * @param ints Integers to store in a buffer.
   * @param bufferN The buffer to store the integers in.
   * @return A buffer containing the specified integers.
   */
  protected IntBuffer getIntBuffer(Collection<Integer> ints, int bufferN) {
    assert (bufferN < 3);
    assert (bufferN >= 0);
    IntBuffer buffer = getBuffer(bufferN, ints.size());
    int index = 0;
    for (Integer i : ints) {
      buffer.put(index, i);
      index++;
    }
    return buffer;
  }

  /**
   * Add a unit clause to the solver.
   *
   * @param a The literal to be asserted true in the solver.
   * @return False if the formula is trivially unsatisfiable after adding this clause (or if it was
   *     already trivially unsatisfiable), else returns true.
   */
  public boolean addClause(Lit a) {
    validate(a);
    return MonosatJNI.addUnitClause(getSolverPtr(), a.l);
  }

  // Solver API

  /**
   * Add a binary clause to the solver.
   *
   * @param a The first literal of the clause.
   * @param b The second literal of the clause.
   * @return False if the formula is trivially unsatisfiable after adding this clause (or if it was
   *     already trivially unsatisfiable), else returns true.
   */
  public boolean addClause(Lit a, Lit b) {
    validate(a);
    validate(b);
    return MonosatJNI.addBinaryClause(getSolverPtr(), a.l, b.l);
  }

  /**
   * Add a ternary clause to the solver.
   *
   * @param a The first literal of the clause.
   * @param b The second literal of the clause.
   * @param c The third literal of the clause.
   * @return False if the formula is trivially unsatisfiable after adding this clause (or if it was
   *     already trivially unsatisfiable), else returns true.
   */
  public boolean addClause(Lit a, Lit b, Lit c) {
    validate(a);
    validate(b);
    validate(c);
    return MonosatJNI.addTertiaryClause(getSolverPtr(), a.l, b.l, c.l);
  }

  /**
   * Add a clause to the solver. Returns false if the formula is trivially unsatisfiable after
   * adding this clause (or if it was already trivially unsatisfiable), else returns true. If args
   * is empty, then the solver will be trivially unsatisfiable.
   *
   * @param args The literals to add as a clause to the solver (forcing at least one of them to be
   *     true)
   * @return False if, after adding the clause, the solver is trivially unsatisfiable; true
   *     otherwise.
   */
  public boolean addClause(Lit... args) {
    validate(args);
    return MonosatJNI.addClause(getSolverPtr(), getLitBuffer(args, 1), args.length);
  }

  /**
   * Add a clause to the solver. Returns false if the formula is trivially unsatisfiable after
   * adding this clause (or if it was already trivially unsatisfiable), else returns true. If args
   * is empty, then the solver will be trivially unsatisfiable.
   *
   * @param args The clause to add to the solver
   * @return False if, after adding the clause, the solver is trivially unsatisfiable; true
   *     otherwise.
   */
  public boolean addClause(Collection<Lit> args) {
    validate(args);
    return MonosatJNI.addClause(getSolverPtr(), getLitBuffer(args, 1), args.size());
  }

  /**
   * Either find a satisfying solution to the constraints in te formula, or prove the formula to be
   * unsatisfiable.
   *
   * @return True if a satisfying solution was found, false otherwise.
   */
  public boolean solve() {
    return MonosatJNI.solve(getSolverPtr());
  }

  // basic this functions

  /**
   * Either find a satisfying solution to the constraints in te formula, or prove the formula to be
   * unsatisfiable, while temporarily enforcing the literals in assumptions to be true.
   *
   * @return True if a satisfying solution was found, false otherwise.
   */
  public boolean solve(Lit... assumptions) {
    validate(assumptions);
    return MonosatJNI.solveAssumptions(getSolverPtr(), getLitBuffer(assumptions, 0), assumptions.length);
  }

  /**
   * Either find a satisfying solution to the constraints in te formula, or prove the formula to be
   * unsatisfiable, while temporarily enforcing the literals in assumptions to be true.
   *
   * @return True if a satisfying solution was found, false otherwise.
   */
  public boolean solve(Collection<Lit> assumptions) {
    validate(assumptions);
    return MonosatJNI.solveAssumptions(getSolverPtr(), getLitBuffer(assumptions), assumptions.size());
  }

  /**
   * Sets the (approximate) time limit in seconds before returning empty from
   * solveLimited(); ignored by solve(). Set to <0 to disable time limit. Note that the solver will
   * only respect this limit on a best effort basis.
   *
   * Note: setTimeLimit has no effect on mac.
   *
   * @param seconds The (approximate) number of seconds to limit the solver to, in the next call to
   *     solveLimited().
   */
  public void setTimeLimit(int seconds) {
    MonosatJNI.setTimeLimit(getSolverPtr(), seconds);
  }

  /**
   * Sets the maximum number of (further) conflicts allowed in the solver before returning
   * empty from solveLimited(); ignored by solve(). Set to <0 to disable conflict limit.
   *
   * @param conflicts The number of additional conflicts rounds allowed by the next call to
   *     solveLimited().
   */
  public void setConflictLimit(int conflicts) {
    MonosatJNI.setConflictLimit(getSolverPtr(), conflicts);
  }

  /**
   * Get the total number of conflicts that have occurred in the solver.
   * @return The total number of conflicts that have occurred in the solver.
   */
  public long nConflicts(){
    return MonosatJNI.nConflicts(getSolverPtr());
  }

  /**
   * Get the total number of unit propagations that have occurred in the solver.
   * @return The total number of unit propagations that have occurred in the solver.
   */
  public long nPropagations(){
    return MonosatJNI.nPropagations(getSolverPtr());
  }

  /**
   * Sets the maximum number of (further) propagation rounds allowed in the solver before returning
   * empty from solveLimited(); ignored by solve(). Set to <0 to disable propagation
   * limit.
   *
   * @param propagations The number of additional propagation rounds allowed by the next call to
   *     solveLimited().
   */
  public void setPropagationLimit(int propagations) {
    MonosatJNI.setPropagationLimit(getSolverPtr(), propagations);
  }

  /**
   * Attempt to find a satisfying solution, but return Optional.empty() if any resource limits are
   * violated. To set resource limits, see:
   *
   * @see #setTimeLimit(int)
   * @see #setConflictLimit(int)
   * @see #setPropagationLimit(int)
   * @return An Optional containing True if a satisfying solution was found, False if the
   *     constraints were proven UNSAT, and Optional.empty() if the solver could neither find a
   *     solution nor prove the constraints UNSAT within the resource limits.
   */
  public Optional<Boolean> solveLimited() {
    int result = MonosatJNI.solveLimited(getSolverPtr());
    assert (result >= 0);
    assert (result <= 2);
    return LBool.toLbool(result).toOpt();
  }

  /**
   * Attempt to find a satisfying solution, but return Optional.empty() if any resource limits are
   * violated. To set resource limits, see:
   *
   * @see #setTimeLimit(int)
   * @see #setConflictLimit(int)
   * @see #setPropagationLimit(int)
   * @param assumptions Literals to enforce temporarily during the solve call.
   * @return An Optional containing True if a satisfying solution was found, False if the
   *     constraints were proven UNSAT, and Optional.empty() if the solver could neither find a
   *     solution nor prove the constraints UNSAT within the resource limits.
   */
  public Optional<Boolean> solveLimited(Collection<Lit> assumptions) {
    validate(assumptions);
    int result =
        MonosatJNI.solveAssumptionsLimited(
            getSolverPtr(), getLitBuffer(assumptions), assumptions.size());
    assert (result >= 0);
    assert (result <= 2);
    return LBool.toLbool(result).toOpt();
  }

  /**
   * Attempt to find a satisfying solution, but return Optional.empty() if any resource limits are
   * violated. To set resource limits, see:
   *
   * @see #setTimeLimit(int)
   * @see #setConflictLimit(int)
   * @see #setPropagationLimit(int)
   * @param assumptions Literals to enforce temporarily during the solve call.
   * @return An Optional containing True if a satisfying solution was found, False if the
   *     constraints were proven UNSAT, and Optional.empty() if the solver could neither find a
   *     solution nor prove the constraints UNSAT within the resource limits.
   */
  public Optional<Boolean> solveLimited(Lit... assumptions) {
    validate(assumptions);
    int result =
        MonosatJNI.solveAssumptionsLimited(
            getSolverPtr(), getLitBuffer(assumptions, 0), assumptions.length);
    assert (result >= 0);
    assert (result <= 2);
    return LBool.toLbool(result).toOpt();
  }

  /**
   * Returns true if the last solve() or solveLimited() call found an optimal solution. This is
   * normally always the case, unless resource limits are enforced.
   *
   * @return True if the most recent solve() or solveLimited() call found an optimal solution.
   */
  public boolean lastSolutionWasOptimal() {
    return MonosatJNI.lastSolutionWasOptimal(getSolverPtr());
  }

  /**
   * If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a
   * subset of the assumptions which are sufficient to cause the instance to be UNSAT).
   *
   * @return A set of literals, selected from among the assumptions in the most recent solve() call,
   *     at least one of which must be true in any satisfying solution.
   */
  public List<Lit> getConflictClause() {
    return getConflictClause(false);
  }

  /**
   * If the last solve call was UNSAT due to one or more assumption literals, this method returns a
   * subset of those assumption literals that are sufficient to keep make the instance UNSAT.
   *
   * @param minimize If true, a (possibly expensive) search will be applied in the solver to find a
   *     locally minimal set of literals that are mutually UNSAT. Else, an inexpensive, best effort
   *     set of literals will be returned.
   * @return A set of literals, selected from among the assumptions in the most recent solve() call,
   *     at least one of which must be true in any satisfying solution.
   */
  public List<Lit> getConflictClause(boolean minimize) {
    ArrayList<Lit> store = new ArrayList<Lit>();
    if (minimize) {
      minimizeConflictClause();
    }

    // If the last solution was unsat, then this get the 'conflict clause' produced by the solver (a
    // subset of the assumptions which are sufficient to cause the instance to be UNSAT).
    // Fills the given pointer with the first max_store_size literals of the conflict clause, and
    // returns the number of literals in the conflict clause. Set store_clause to null and
    // max_store_size to 0 to find the size of the conflict clause
    // Returns -1 if the solver has no conflict clause from the most recent solve() call (because
    // that call was not UNSAT)
    int conflict_size = MonosatJNI.getConflictClause(getSolverPtr(), null, 0);
    IntBuffer buf = getBuffer(0, conflict_size);
    int sz = MonosatJNI.getConflictClause(getSolverPtr(), buf, conflict_size);
    assert (sz == conflict_size);
    for (int i = 0; i < conflict_size; i++) {
      int l = buf.get(i);
      Lit lit = toLit(l);
      store.add(lit);
    }
    validate(store);
    return store;
  }

  /**
   * Internal method used to cache newly instantiated literals, so that literals returned from C++
   * can reuse the same objects.
   *
   * @param l The literal to add to the solver.
   */
  protected void registerLit(Lit l) {
    registerLit(l, null);
  }

  /**
   * Internal method used to maintain an internal list of bitvectors in the solver.
   *
   * @param bv The BitVector to add to the solver.
   */
  protected void registerBitVector(BitVector bv) {
    int id = bv.id;
    while (id >= bvmap.size()) {
      bvmap.add(null);
    }
    if (bvmap.get(id) == null) {
      bvmap.set(id, bv);
      allBVs.add(bv);
    }
  }
  /**
   * Internal method to convert C++ bitvector integers into java BitVectors. Will return the
   * existing BitVector object for this ID if one already exists, else, will create a new BitVector
   * object for this id.
   *
   * @param bvID A bitvector id to convert into a BitVector
   * @return A BitVectors object representing the bv.
   */
  protected BitVector toBitVector(int bvID) {
    assert (bvID >= 0);

    assert (bvID
        < nBitVectors()); // the bitvector must have already been declared in the sat solver before
    // this

    while (bvID >= bvmap.size()) {
      bvmap.add(null);
    }

    if (bvmap.get(bvID) == null) {
      BitVector bv = new BitVector(this, this, bvID);
      assert (bvmap.get(bvID) == bv);
      return bv;
    } else {
      return bvmap.get(bvID);
    }
  }
  /**
   * Internal method used to cache newly instantiated literals, so that literals returned from C++
   * can reuse the same objects.
   *
   * @param l The literal to add to the solver.
   */
  protected void registerLit(Lit l, Lit notL) {
    assert (l.l >= 0);
    int literal = l.toInt();
    assert (literal >= 0);
    int var = literal / 2;
    if (notL == null) {
      assert (var
          < nVars()); // the variable must have already been declared in the sat solver before this
    }
    // call
    while (var * 2 + 1 >= allLits.size()) {
      allLits.add(null);
    }

    assert (allLits.get(var * 2) == null);
    assert (allLits.get(var * 2 + 1) == null);

    if (notL == null) {
      if (l.sign()) {
        notL = new Lit(this, var * 2);
      } else {
        notL = new Lit(this, var * 2 + 1);
      }
    }

    if (l.sign()) {
      assert (notL.l == l.l - 1);
      allLits.set(var * 2, notL);
      allLits.set(var * 2 + 1, l);
    } else {
      assert (notL.l == l.l + 1);
      allLits.set(var * 2, l);
      allLits.set(var * 2 + 1, notL);
    }

    assert (!allLits.get(var * 2).sign());
    assert (allLits.get(var * 2 + 1).sign());

    assert (allLits.get(l.l) != null);
    assert (allLits.get(l.l) == l);

    if (l.sign()) {
      positiveLiterals.add(notL);
    } else {
      positiveLiterals.add(l);
    }
  }

  /**
   * Internal method to convert C++ literal integers into java Lits
   *
   * @param literal A literal to convert into a Lit
   * @return A Lit object representing the literal.
   */
  protected Lit toLit(int literal) {
    assert (literal >= 0);
    int var = literal / 2;
    assert (var
        < nVars()); // the variable must have already been declared in the sat solver before this
    // call
    while (var * 2 + 1 >= allLits.size()) {
      allLits.add(null);
    }
    /* if (!MonosatJNI.hasVariable(var)) {
        throw new RuntimeException("No such variable in solver: " + var);
    }*/
    if (allLits.get(var * 2) == null) {
      assert (allLits.get(var * 2 + 1) == null);
      Lit l = new Lit(this, var * 2);
      Lit notL = new Lit(this, var * 2 + 1);
      allLits.set(var * 2, l);
      allLits.set(var * 2 + 1, notL);
      assert (!allLits.get(var * 2).sign());
      assert (allLits.get(var * 2 + 1).sign());

      positiveLiterals.add(l);
    }
    assert (allLits.get(literal) != null);
    assert (allLits.get(literal).l == literal);
    return allLits.get(literal);
  }

  /**
   * Get a collection of all the literals in this solver.
   *
   * @return unmodifiable view of the positive literals in the solver.
   */
  public Collection<Lit> getLits() {
    return positiveLiterals;
  }

  /**
   * Get a collection of all the BitVectors in this solver.
   *
   * @return unmodifiable view of the BitVectors in the solver.
   */
  public Collection<BitVector> getBitVectors() {
    return Collections.unmodifiableList(allBVs);
  }

  /**
   * Given a set of assumptions which are mutually UNSAT, find a locally minimal subset that remains
   * UNSAT. (leaves the original set intact if the literals are not mutually UNSAT)
   *
   * @param literals A set of mutually unsatisfiable literals to minimize.
   * @return A subset of literals that is still UNSAT.
   */
  public List<Lit> minimizeUnsatCore(Lit... literals) {
    validate(literals);
    return minimizeUnsatCore(Arrays.asList(literals));
  }

  // Optimization API

  /**
   * Given a set of assumptions which are mutually UNSAT, find a locally minimal subset that remains
   * UNSAT. (leaves the original set intact if the literals are not mutually UNSAT)
   *
   * @param literals A set of mutually unsatisfiable literals to minimize.
   * @return A subset of literals that is still UNSAT.
   */
  public List<Lit> minimizeUnsatCore(Collection<Lit> literals) {
    validate(literals);
    ArrayList<Lit> store = new ArrayList<Lit>();
    IntBuffer buf = getLitBuffer(literals);
    int core_size = MonosatJNI.minimizeUnsatCore(getSolverPtr(), buf, literals.size());
    assert (core_size >= 0);
    assert (core_size <= literals.size());
    for (int i = 0; i < core_size; i++) {
      int l = buf.get(i);
      Lit lit = toLit(l);
      store.add(lit);
    }
    validate(store);
    return store;
  }

  /**
   * After UNSAT solve calls with assumptions, the solver will find a 'conflict clause' consisting
   * of a subset of the assumptions //which are sufficient to make the solver UNSAT (see
   * getConflictClause). Normally, the conflict clause is produced as a side effect of proving the
   * query unsat, with the solver only removing literals from the conflict clause on a best effort
   * basis. This method will make repeated (and potentially expensive) calls to the SAT solver to
   * attempt to remove further literals from the conflict clause. Afterward, the conflict clause can
   * be obtained using getConflictClause(). NOTE: this function may be expensive, is not required to
   * get a conflict clause; getConflictClause() can be used after any UNSAT call with assumptions,
   * even without calling minimizeConflictClause(). Note also that if any of the setXXXXXLimits()
   * are applied to the solver, this may not produce a locally minimal conflict clause.
   */
  private void minimizeConflictClause() {
    MonosatJNI.minimizeConflictClause(getSolverPtr());
  }

  /** Clear any optimization objectives in the solver. */
  public void clearOptimizationObjectives() {
    MonosatJNI.clearOptimizationObjectives(getSolverPtr());
  }

  /**
   * Add an optimization objective to the solver: maximize the specified BitVector.
   *
   * @param bv The BitVector to maximize in the next solve() or solveLimited() call.
   */
  public void maximizeBV(BitVector bv) {
    validate(bv);
    MonosatJNI.maximizeBV(getSolverPtr(), bvPtr, bv.id);
  }

  /**
   * Add an optimization objective to the solver: minimize the specified BitVector.
   *
   * @param bv The BitVector to minimize in the next solve() or solveLimited() call.
   */
  public void minimizeBV(BitVector bv) {
    validate(bv);
    MonosatJNI.minimizeBV(getSolverPtr(), bvPtr, bv.id);
  }

  /**
   * Add an optimization objective to the solver: maximize the number of true literals from among
   * the given literals.
   *
   * @param literals A collection of literals to maximize in the next solve() or solveLimited()
   *     call.
   */
  public void maximizeLits(Collection<Lit> literals) {
    validate(literals);
    MonosatJNI.maximizeLits(getSolverPtr(), getLitBuffer(literals), literals.size());
  }

  /**
   * Add an optimization objective to the solver: minimize the number of true literals from among
   * the given literals.
   *
   * @param literals A collection of literals to minimize in the next solve() or solveLimited()
   *     call.
   */
  public void minimizeLits(Collection<Lit> literals) {
    validate(literals);
    MonosatJNI.minimizeLits(getSolverPtr(), getLitBuffer(literals), literals.size());
  }

  /**
   * Add an optimization objective to the solver: maximize the weighted number of true literals from
   * among the given literals.
   *
   * @param literals A collection of weighted literals to maximize in the next solve() or
   *     solveLimited() call.
   * @param weights Weights for each literal in literals.
   * @throws IllegalArgumentException If weights.size() != literals.size();
   */
  public void maximizeWeightedLits(List<Lit> literals, List<Integer> weights) {
    validate(literals);
    if (literals.size() != weights.size()) {
      throw new IllegalArgumentException("literals and weights must be the same length");
    }
    MonosatJNI.maximizeWeightedLits(
        getSolverPtr(), getLitBuffer(literals), getIntBuffer(weights, 1), literals.size());
  }

  /**
   * Add an optimization objective to the solver: minimize the weighted number of true literals from
   * among the given literals.
   *
   * @param literals A collection of weighted literals to minimize in the next solve() or
   *     solveLimited() call.
   * @param weights Weights for each literal in literals.
   * @throws IllegalArgumentException If weights.size() != literals.size();
   */
  public void minimizeWeightedLits(List<Lit> literals, List<Integer> weights) {
    validate(literals);
    if (literals.size() != weights.size()) {
      throw new IllegalArgumentException("literals and weights must be the same length");
    }
    MonosatJNI.minimizeWeightedLits(
        getSolverPtr(), getLitBuffer(literals), getIntBuffer(weights, 1), literals.size());
  }

  /**
   * Enforce that at most one of the specified literals may be true. If 1 or fewer arguments are
   * given, has no effect. If exactly 2 arguments are given, this is the same as:
   * assertOr(args[0],args[1])
   *
   * @param args An array of Lits, at most one of which can be true.
   */
  public void assertAtMostOne(Lit... args) {
    validate(args);
    assertAtMostOne(Arrays.asList(args));
  }

  /**
   * Enforce that at most one of the specified literals may be true. If 1 or fewer arguments are
   * given, has no effect. If exactly 2 arguments are given, this is the same as:
   * assertOr(args[0],args[1])
   *
   * @param args A collection of Lits, at most one of which can be true.
   */
  public void assertAtMostOne(Collection<Lit> args) {
    // simple at-most-one constraint: asserts that at most one of the set of literals may be true.
    validate(args);
    if (args.size() <= 6) {
      // make this constant configurable in the future
      // for small enough sets of literals, directly instantiate the binary clauses constraining
      // them rather than introduce an amo theory
      MonosatJNI.AssertAMO(getSolverPtr(), getLitBuffer(args), args.size());
    } else {
      MonosatJNI.at_most_one_lit(getSolverPtr(), getLitBuffer(args, 0), args.size());
    }
  }

  /**
   * Enforce a pseudo-Boolean constraint. The number of true literals from among args must satisfy
   * comparison c relative to compare to.
   *
   * <p>For example, if c is Comparison.LEQ, and compareTo is 3, then at most 3 literals from args
   * may be true.
   *
   * @param args One or more literals, whose sum will be compared.
   * @param c The comparison operation to perform.
   * @param compareTo The constant to compare the sum of the true lits in args to.
   */
  public void assertPB(Collection<Lit> args, Comparison c, int compareTo) {
    assertPB(args, null, c, compareTo);
  }

  /**
   * Enforce a weighted pseudo-Boolean constraint. The weighted number of true literals from among
   * args must satisfy comparison c relative to compare to.
   *
   * @param args One or more literals, whose sum will be compared.
   * @param weights Weights for each literal in args (if fewer weights than args are supplied, the
   *     remaining weights will be set to '1').
   * @param c The comparison operation to perform.
   * @param compareTo The constant to compare the sum of the true lits in args to.
   */
  public void assertPB(Collection<Lit> args, Collection<Integer> weights, Comparison c, int compareTo) {
    ArrayList<Lit> a = new ArrayList<Lit>();
    ArrayList<Integer> w = null;
    args.forEach(a::add);
    if(weights!=null){
      w = new ArrayList<Integer>();
      weights.forEach(w::add);
    }
    assertPB(a,w, c, compareTo);
  }

  /**
   * Enforce a weighted pseudo-Boolean constraint. The weighted number of true literals from among
   * args must satisfy comparison c relative to compare to.
   *
   * @param args A collection of args, whose sum will be compared.
   * @param weights Weights for each literal in args (if fewer weights than args are supplied, the
   *     remaining weights will be set to '1').
   * @param c The comparison operation to perform.
   * @param compareTo The constant to compare the sum of the true lits in args to.
   */
  public void assertPB(List<Lit> args, List<Integer> weights, Comparison c, int compareTo) {
    validate(args);


    IntBuffer wt_buffer = getBuffer(1, args.size());
    int n_wts = 0;
    if (weights != null) {
      if (weights.size() > args.size()) {
        throw new IllegalArgumentException("Must not have more weights then literals.");
      }
      n_wts = Math.min(args.size(), weights.size());
      int i = 0;
      for (Integer w : weights) {
        wt_buffer.put(i, w);
        if (++i >= weights.size()) {
          break;
        }
      }
    }
    for (int i = n_wts; i < args.size(); i++) {
      wt_buffer.put(i, 1); // default weight is 1
    }
    switch (c) {
      case LT:
        MonosatJNI.assertPB_lt(getSolverPtr(), compareTo, args.size(), getLitBuffer(args), wt_buffer);
        break;
      case LEQ:
        MonosatJNI.assertPB_leq(getSolverPtr(), compareTo, args.size(), getLitBuffer(args), wt_buffer);
        break;
      case EQ:
        MonosatJNI.assertPB_eq(getSolverPtr(), compareTo, args.size(), getLitBuffer(args), wt_buffer);
        break;
      case GEQ:
        MonosatJNI.assertPB_geq(getSolverPtr(), compareTo, args.size(), getLitBuffer(args), wt_buffer);
        break;
      case GT:
        MonosatJNI.assertPB_gt(getSolverPtr(), compareTo, args.size(), getLitBuffer(args), wt_buffer);
        break;
      case NEQ:
      default:
        throw new UnsupportedOperationException("Unsupported comparison " + c);
    }
  }

  /**
   * Immediately convert any pseudo Boolean constraints in the solver into clauses. This function
   * does not need to be manually called, as it will be called automatically before 'solve' calls.
   */
  protected void flushPB() {
    MonosatJNI.flushPB(getSolverPtr());
  }

  /**
   * Create a constant-valued BitVector of the given width. Caches constant bitvectors of value <=
   * MAX_CACHE_CONST
   *
   * @param width The width of the BitVector to create.
   * @param constant The constant value of the BitVector.
   * @return A new BitVector of the corresponding size.
   */
  public BitVector bv(int width, long constant) {
    if (constant >= 0 && constant <= MAX_CACHE_CONST) {
      while (cached_bvs.size() <= width) {
        cached_bvs.add(new ArrayList<>());
      }

      int small_const = (int) constant;
      while (cached_bvs.get(width).size() <= constant) {
        cached_bvs.get(width).add(null);
      }
      if (cached_bvs.get(width).get(small_const) == null) {
        cached_bvs.get(width).set(small_const, new BitVector(this, width, constant));
      }
      assert (cached_bvs.get(width).get(small_const) != null);
      return cached_bvs.get(width).get(small_const);
    }
    return new BitVector(this, width, constant);
  }

    /**
     * Assign a name to a literal.
     * <p>
     * Literals in MonoSAT may have zero or more names.
     * Positive and negative signed literals of a variable share the same names.
     * @param l The literal to add a name to.
     * @param name The name to assign to the literal.
     * If the string is empty, then no name will be
     * added to this literal. Otherwise, the name must be unique,
     * and is restricted to printable ASCII characters.
     *
     * @throws IllegalArgumentException If the new name is not a valid ID.
     */
    public void addName(Lit l, String name) {
        validate(l);
        if(name.contains("~")){
            throw new IllegalArgumentException("Literal IDs may not include \"~\":" + name + "\" is not a valid name.");
        }
        if(name!=null && name.length()>0) {
            String proposedName = MonosatJNI.validID(name);
            MonosatJNI.addLiteralName(getSolverPtr(), l.l, proposedName);
            if (l._name == null || l._name.length() == 0) {
                l._name = proposedName;
            }
            if (l.not()._name == null || l.not()._name.length() == 0) {
                l.not()._name = "~" + proposedName;
            }
        }
    }
    /**
     * Assign a name to a BitVector.
     * <p>
     * BitVectors in MonoSAT may have zero or more names.
     * @param bv The BitVector to add a name to.
     * @param name The name to assign to the BitVector.
     * If the string is empty, then no name will be
     * added to this BitVector. Otherwise, the name must be unique,
     * and is restricted to printable ASCII characters.
     *
     * @throws IllegalArgumentException If the new name is not a valid ID.
     */
    public void addName(BitVector bv, String name) {
        validate(bv);
        if(name!=null && name.length()>0) {
            String proposedName = MonosatJNI.validID(name);
            MonosatJNI.setBitvectorName(getSolverPtr(),bvPtr, bv.getID(), proposedName);
            if (bv._name == null || bv._name.length() == 0) {
                bv._name = proposedName;
            }
        }
    }
  /**
   * Returns an iterator over the (positive) literals in the solver. Each literal in the iterator
   * has positive sign.
   *
   * @return An iterator over the literals in the solver.
   */
  public Iterable<Lit> literals() {
    return new LitIterator(false);
  }

  /**
   * Returns an iterator over the named literals in the solver. Each literal in the iterator has
   * positive sign.
   *
   * @return An iterator over the named literals in the solver.
   */
  public Iterable<Lit> namedLiterals() {
    return new LitIterator(true);
  }

  /** An iterator over literals in the solver. */
  public class LitIterator implements java.util.Iterator<Lit>,java.lang.Iterable<Lit> {
    private int index = 0;

    final boolean named;

    protected LitIterator(boolean named) {
      this.named = named;
    }

    @Override
    public boolean hasNext() {
      if (named) {
        return index < MonosatJNI.nNamedLiterals(Solver.this.getSolverPtr());
      } else {
        return index < Solver.this.nVars();
      }
    }

    @Override
    public Lit next() {
      if (!hasNext()) {
        throw new IndexOutOfBoundsException();
      }
      if (named) {
        return toLit(MonosatJNI.getNamedLiteralN(Solver.this.getSolverPtr(), index++));
      } else {
        return toLit((index++) * 2);
      }
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<Lit> iterator() {
      return this;
    }
  }

  /**
   * Returns an iterator over the (positive) literals in the solver. Each literal in the iterator
   * has positive sign.
   *
   * @return An iterator over the literals in the solver.
   */
  public Iterable<BitVector> bitvectors() {
    return new BVIterator(false);
  }

  /**
   * Returns an iterator over the named literals in the solver. Each literal in the iterator has
   * positive sign.
   *
   * @return An iterator over the named literals in the solver.
   */
  public Iterable<BitVector> namedBitVectors() {
    return new BVIterator(true);
  }

  /** An iterator over bitvectors in the solver. */
  public class BVIterator implements java.util.Iterator<BitVector>,java.lang.Iterable<BitVector> {
    private int index = 0;

    final boolean named;

    protected BVIterator(boolean named) {
      this.named = named;
    }

    @Override
    public boolean hasNext() {
      if (named) {
        return index < MonosatJNI.nNamedBitvectors(Solver.this.getSolverPtr(), Solver.this.bvPtr);
      } else {
        return index < Solver.this.nBitVectors();
      }
    }

    @Override
    public BitVector next() {
      if (!hasNext()) {
        throw new IndexOutOfBoundsException();
      }
      if (named) {
        return toBitVector(
            MonosatJNI.getNamedBitvectorN(Solver.this.getSolverPtr(), Solver.this.bvPtr, index++));
      } else {
        return toBitVector(index++);
      }
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<BitVector> iterator() {
      return this;
    }
  }
  /**
   * Retrieve an existing named literal from the solver, by looking up its name. Note that the
   * returned literal will never be negated.
   *
   * @param name The name of the literal to check for.
   * @return The matching literal, if it exists.
   * @throws IllegalArgumentException If there is no literal in the solver with this name (or if
   *     name is empty).
   */
  public Lit getLiteral(String name) {
    if (name == null || name.length() == 0) {
      throw new IllegalArgumentException("Name must not be empty");
    }
    if (name == "True") {
      return Lit.True;
    } else if (name == "False") {
      return Lit.False;
    }

    int literal = MonosatJNI.getLiteral(getSolverPtr(), MonosatJNI.validID(name));

    if (literal >= 0) {
      Lit lit = toLit(literal);
      validate(lit);
      return lit;
    } else {
      throw new IllegalArgumentException("No variable with name " + name);
    }
  }

  /**
 * Test if the solver has an existing named literal with the given name.
 * Note that an IllegalArgumentException will be thrown if the name is not a valid
 * identifier.
 *
 * @param name The name of the literal to check for.
 * @return True if there is a literal with the given name in the solver.
 * @throws IllegalArgumentException If the name is not a valid identifier.
 */
  public boolean hasLiteral(String name) {
      if (name == null || name.length() == 0) {
          return false;
      }
      return MonosatJNI.hasLiteralWithName(getSolverPtr(), MonosatJNI.validID(name));
  }

  /**
   * Retrieve an existing named bitvector from the solver, by looking up its name.
   *
   * @param name The name of the bitvector to check for.
   * @return The matching bitvector, if it exists.
   * @throws IllegalArgumentException If there is no bitvector in the solver with this name (or if
   *     name is empty).
   */
  public BitVector getBitVector(String name) {
    if (name == null || name.length() == 0) {
      throw new IllegalArgumentException("Name must not be empty");
    }
    int bvID = MonosatJNI.getBitvector(getSolverPtr(), bvPtr, MonosatJNI.validID(name));

    if (bvID >= 0) {
      return new BitVector(this, this, bvID);
    } else {
      throw new IllegalArgumentException("No Btivector with name " + name);
    }
  }

  /**
   * Retrieve an existing named graph from the solver, by looking up its name.
   *
   * @param name The name of the literal to check for.
   * @return The matching literal, if it exists.
   * @throws IllegalArgumentException If there is no literal in the solver with this name (or if
   *     name is empty).
   */
  public Graph getGraph(String name) {
    if (name == null || name.length() == 0) {
      throw new IllegalArgumentException("Name must not be empty");
    }
    if (name.length() == 0) {
      throw new IllegalArgumentException("Graph string name must have length > 0");
    }
    long graphPtr = MonosatJNI.getGraph(getSolverPtr(), MonosatJNI.validID(name));
    if (graphPtr == 0) {
      throw new IllegalArgumentException("No graph with name " + name);
    }
    if (allGraphs.containsKey(graphPtr)) {
      return allGraphs.get(graphPtr);
    } else {
      Graph g = new Graph(this, graphPtr);
      assert (g.name().equals(name));
      return g;
    }
  }

  /**
   * Retrieve an existing literal from the solver, by looking up its internal integer ID.
   *
   * @param literal The integer of the literal to check for.
   * @return The matching literal, if it exists.
   * @throws IllegalArgumentException If there is no literal in the solver with this integer.
   */
  protected Lit getLiteral(int literal) {
    if (literal < 0) {
      throw new IllegalArgumentException("Invalid literal: " + literal);
    }
    Lit lit = toLit(literal);
    validate(lit);
    return lit;
  }

  /**
   * Reset any decisions in the solver, while preserving the current constraints. There is normally
   * no need to manually call restart(), as this is called automatically before calling solve().
   */
  public void restart() {
    MonosatJNI.backtrack(getSolverPtr());
  }

  /**
   * Check if the solver has a satisfying model to its constraints.
   *
   * @return True if the solver has a satisfying model to its constraints.
   */
  public boolean hasModel() {
    return MonosatJNI.hasModel(getSolverPtr());
  }

  /**
   * Return true if the solver has not yet proven its constraints to be UNSAT. If this returns
   * false, then all future calls to 'solve' will return false.
   *
   * @return true if the solver has not yet proven its constraints to be UNSAT.
   */
  public boolean ok() {
    return MonosatJNI.ok(getSolverPtr());
  }

  /**
   * Return the negation of a literal.
   * <p>
   * Implementation note:<br>
   * Lit.not() is well optimized and inexpensive, in both Java and the native library.
   * Double negations will be eliminated and are essentially free.
   * No new objects are instantiated by Lit.not().
   * It is always the case that <code>a.not().not() == a<code/>
   * @param a The literal to negate.
   * @return The negation of a.
   */
  public Lit not(Lit a) {
    validate(a);
    int l = a.toInt();
    l = l ^ 1; // bit twiddle odd to even
    assert (l >= 0);
    assert (l < allLits.size());
    return allLits.get(l);
  }

  /**
   * Create a new literal that evaluates to true if a is false, or if a is true and b is true.
   *
   * @param a The pre-condition.
   * @param b The post-condition.
   * @return A new literal representing an implication gate.
   */
  public Lit implies(Lit a, Lit b) {
    validate(a, b);
    return or(a.not(), b);
  }

  /**
   * Create a new literal, that is true if a and b both evaluate to the same truth value.
   *
   * @param a The first literal to test equality.
   * @param b The second literal to test equality.
   * @return A literal that is true if a and b both evaluate to the same truth value.
   */
  public Lit equal(Lit a, Lit b) {
    validate(a, b);
    return xnor(a, b);
  }

  /**
   * Create a new literal that will evaluate to 'then' if condition is true, and to 'els' otherwise.
   *
   * @param condition The condition literal to test.
   * @param then The value of the returned literal if 'condition' is true.
   * @param els The value of the returned literal if 'condition' is false.
   * @return A new literal, equal to 'then' if condition is true, and equal to 'els' if condition is
   *     false.
   */
  public Lit ite(Lit condition, Lit then, Lit els) {
    validate(condition, then, els);
    return this.toLit(MonosatJNI.Ite(this.getSolverPtr(), condition.l, then.l, els.l));
  }

  /**
   * Create a new literal that is true if all the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the And gate.
   * @return A new literal, true if all the literals in args are true.
   */
  public Lit and(Lit... args) {
    validate(args);
    if (args.length == 0) {
      return Lit.True;
    } else if (args.length == 1) {
      return args[0];
    } else if (args.length == 2) {

      return this.toLit(MonosatJNI.And(this.getSolverPtr(), args[0].l, args[1].l));
    }
    return and(Arrays.asList(args));
  }

  // Higher level constructs for the solver
  // Literal level constructs

  /**
   * Create a new literal that is true if all the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the And gate.
   * @return A new literal, true if all the literals in args are true.
   */
  public Lit and(Collection<Lit> args) {
    validate(args);
    return this.toLit(MonosatJNI.Ands(this.getSolverPtr(), this.getLitBuffer(args), args.size()));
  }

  /**
   * Create a new literal that is true if any of the literals in args are true. If args is empty,
   * return Lit.False.
   *
   * @param args The arguments to the Or gate.
   * @return A new literal, true if any of the literals in args are true.
   */
  public Lit or(Lit... args) {
    validate(args);
    if (args.length == 0) {
      return Lit.False;
    } else if (args.length == 1) {
      return args[0];
    } else if (args.length == 2) {
      return this.toLit(MonosatJNI.Or(this.getSolverPtr(), args[0].l, args[1].l));
    }
    return or(Arrays.asList(args));
  }

  /**
   * Create a new literal that is true if any of the literals in args are true. If args is empty,
   * return Lit.False.
   *
   * @param args The arguments to the Or gate.
   * @return A new literal, true if any of the literals in args are true.
   */
  public Lit or(Collection<Lit> args) {
    validate(args);
    return this.toLit(MonosatJNI.Ors(this.getSolverPtr(), this.getLitBuffer(args), args.size()));
  }

  /**
   * Create a new literal that is true if not all of the literals in args are true. If args is
   * empty, return Lit.False.
   *
   * @param args The arguments to the Nand gate.
   * @return A new literal, true if not all the literals in args are true.
   */
  public Lit nand(Lit... args) {
    validate(args);
    if (args.length == 0) {
      return this.toLit(MonosatJNI.Nands(this.getSolverPtr(), getBuffer(0, 0), 0));
    } else if (args.length == 1) {
      return args[0].not();
    } else if (args.length == 2) {

      return this.toLit(MonosatJNI.Nand(this.getSolverPtr(), args[0].l, args[1].l));
    }
    return nand(Arrays.asList(args));
  }

  /**
   * Create a new literal that is true if not all of the literals in args are true. If args is
   * empty, return Lit.False.
   *
   * @param args The arguments to the Nand gate.
   * @return A new literal, true if not all the literals in args are true.
   */
  public Lit nand(Collection<Lit> args) {
    validate(args);
    return this.toLit(MonosatJNI.Nands(this.getSolverPtr(), this.getLitBuffer(args), args.size()));
  }

  /**
   * Create a new literal that is true if none of the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the Nor gate.
   * @return A new literal, true if none of the literals in args are true.
   */
  public Lit nor(Lit... args) {
    validate(args);
    if (args.length == 0) {
      return this.toLit(MonosatJNI.Nors(this.getSolverPtr(), getBuffer(0, 0), 0));
    } else if (args.length == 1) {
      return args[0];
    } else if (args.length == 2) {
      return this.toLit(MonosatJNI.Nor(this.getSolverPtr(), args[0].l, args[1].l));
    }
    return nor(Arrays.asList(args));
  }

  /**
   * Create a new literal that is true if none of the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the Nor gate.
   * @return A new literal, true if none of the literals in args are true.
   */
  public Lit nor(Collection<Lit> args) {
    validate(args);
    return this.toLit(MonosatJNI.Nors(this.getSolverPtr(), this.getLitBuffer(args), args.size()));
  }

  /**
   * Create a new literal that is true if an odd number of the literals in args are true. If args is
   * empty, return Lit.True.
   *
   * @param args The arguments to the Xor gate.
   * @return A new literal, true if an odd number of the literals in args are true.
   */
  public Lit xor(Lit... args) {
    validate(args);
    if (args.length == 0) {
      return this.toLit(MonosatJNI.Xors(this.getSolverPtr(), getBuffer(0, 0), 0));
    } else if (args.length == 1) {
      return args[0].not();
    } else if (args.length == 2) {

      return this.toLit(MonosatJNI.Xor(this.getSolverPtr(), args[0].l, args[1].l));
    }
    return xor(Arrays.asList(args));
  }

  /**
   * Create a new literal that is true if an odd number of the literals in args are true. If args is
   * empty, return Lit.True.
   *
   * @param args The arguments to the Xor gate.
   * @return A new literal, true if an odd number of the literals in args are true.
   */
  public Lit xor(Collection<Lit> args) {
    validate(args);
    return this.toLit(MonosatJNI.Xors(this.getSolverPtr(), this.getLitBuffer(args), args.size()));
  }

  /**
   * Create a new literal that is true if an even number of the literals in args are true. If args
   * is empty, return Lit.True.
   *
   * @param args The arguments to the Xnor gate.
   * @return A new literal, true if an even number of the literals in args are true.
   */
  public Lit xnor(Lit... args) {
    validate(args);
    if (args.length == 0) {
      return this.toLit(MonosatJNI.Xnors(this.getSolverPtr(), getBuffer(0, 0), 0));
    } else if (args.length == 1) {
      return args[0];
    } else if (args.length == 2) {

      return this.toLit(MonosatJNI.Xnor(this.getSolverPtr(), args[0].l, args[1].l));
    }
    return xnor(Arrays.asList(args));
  }

  /**
   * Create a new literal that is true if an even number of the literals in args are true. If args
   * is empty, return Lit.True.
   *
   * @param args The arguments to the Xnor gate.
   * @return A new literal, true if an even number of the literals in args are true.
   */
  public Lit xnor(Collection<Lit> args) {
    validate(args);
    return this.toLit(MonosatJNI.Xnors(this.getSolverPtr(), this.getLitBuffer(args), args.size()));
  }

  /**
   * Create a new literal that evaluates to true if a is false, or if a is true and if at least one element of args
   * it true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public Lit impliesOr(Lit a, Collection<Lit> args) {
    validate(a);
    validate(args);
    return this.toLit(MonosatJNI.ImpliesOr_(this.getSolverPtr(), this.getLitBuffer(args), args.size(),a.toInt()));
  }

  /**
   * Create a new literal that evaluates to true if a is false, or if a is true and if at least one element of args
   * it true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public Lit impliesOr(Lit a, Lit... args) {
    validate(a);
    validate(args);
    return this.toLit(MonosatJNI.ImpliesOr_(this.getSolverPtr(), this.getLitBuffer(args,0), args.length,a.toInt()));
  }

  /**
   * Create a new literal that evaluates to true if a is false, or if a is true and if at least one element of args
   * it true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public Lit impliesAnd(Lit a, Collection<Lit> args) {
    validate(a);
    validate(args);
    return this.toLit(MonosatJNI.ImpliesAnd(this.getSolverPtr(), this.getLitBuffer(args), args.size(),a.toInt()));
  }

  /**
   * Create a new literal that evaluates to true if a is false, or if a is true and if at least one element of args
   * it true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public Lit impliesAnd(Lit a, Lit... args) {
    validate(a);
    validate(args);
    return this.toLit(MonosatJNI.ImpliesAnd(this.getSolverPtr(), this.getLitBuffer(args,0), args.length,a.toInt()));
  }
  /**
   * Assert that Lit a must be true in the solver.
   *
   * @param a The literal to assert true.
   */
  public void assertTrue(Lit a) {
    validate(a);
    MonosatJNI.Assert(this.getSolverPtr(), a.l);
  }

  // Assertion forms.

  /**
   * Assert that Lit a must be false in the solver.
   *
   * @param a The literal to assert false.
   */
  public void assertFalse(Lit a) {
    validate(a);
    MonosatJNI.Assert(this.getSolverPtr(), a.not().l);
  }

  /**
   * Assert that all of the literals in args are true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert true.
   */
  public void assertAnd(Lit... args) {
    validate(args);
    if (args.length == 0) {
      // do nothing
    } else if (args.length == 1) {
      assertTrue(args[0]);
    } else if (args.length == 2) {
      MonosatJNI.AssertAnd(this.getSolverPtr(), args[0].l, args[1].l);
    } else {
      assertAnd(Arrays.asList(args));
    }
  }

  /**
   * Assert that all of the literals in args are true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert true.
   */
  public void assertAnd(Collection<Lit> args) {
    validate(args);
    MonosatJNI.AssertAnds(this.getSolverPtr(), this.getLitBuffer(args), args.size());
  }

  /**
   * Assert that at least one of the literals in args are true. Trivially contradictory if args is
   * empty.
   *
   * <p>Implementation note: This call is equivalent to Solver.addClause().
   *
   * @param args The arguments to assert at least one must hold.
   */
  public void assertOr(Lit... args) {
    validate(args);
    if (args.length == 0) {
      // this is a contradiction
      MonosatJNI.AssertOrs(this.getSolverPtr(), getBuffer(0, 0), 0);
    } else if (args.length == 1) {
      assertTrue(args[0]);
    } else if (args.length == 2) {
      MonosatJNI.AssertOr(this.getSolverPtr(), args[0].l, args[1].l);
    } else {
      assertOr(Arrays.asList(args));
    }
  }

  /**
   * Assert that at least one of the literals in args are true. Trivially contradictory if args is
   * empty.
   *
   * <p>Implementation note: This call is equivalent to Solver.addClause().
   *
   * @param args The arguments to assert at least one must hold.
   */
  public void assertOr(Collection<Lit> args) {
    validate(args);
    MonosatJNI.AssertOrs(this.getSolverPtr(), this.getLitBuffer(args), args.size());
  }

  /**
   * Assert that at least one of the given args must be false. Trivially contradictory if args is
   * empty.
   *
   * @param args The arguments to assert at least one must be false.
   */
  public void assertNand(Lit... args) {
    validate(args);
    if (args.length == 0) {
      MonosatJNI.AssertNands(this.getSolverPtr(), getBuffer(0, 0), 0);
    } else if (args.length == 1) {
      assertTrue(args[0].not());
    } else if (args.length == 2) {
      MonosatJNI.AssertNand(this.getSolverPtr(), args[0].l, args[1].l);
    } else {
      assertNand(Arrays.asList(args));
    }
  }

  /**
   * Assert that at least one of the given args must be false. Trivially contradictory if args is
   * empty.
   *
   * @param args The arguments to assert at least one must be false.
   */
  public void assertNand(Collection<Lit> args) {
    validate(args);
    MonosatJNI.AssertNands(this.getSolverPtr(), this.getLitBuffer(args), args.size());
  }

  /**
   * Assert that none of the given args may be true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert false.
   */
  public void assertNor(Lit... args) {
    validate(args);
    if (args.length == 0) {
      MonosatJNI.AssertNors(this.getSolverPtr(), getBuffer(0, 0), 0);
    } else if (args.length == 1) {
      assertTrue(args[0].not());
    } else if (args.length == 2) {
      MonosatJNI.AssertNor(this.getSolverPtr(), args[0].l, args[1].l);
    } else {
      assertNor(Arrays.asList(args));
    }
  }

  /**
   * Assert that none of the given args may be true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert false.
   */
  public void assertNor(Collection<Lit> args) {
    validate(args);
    MonosatJNI.AssertNors(this.getSolverPtr(), this.getLitBuffer(args), args.size());
  }

  /**
   * Assert that an odd number of args must hold. Trivially contradictory if args is empty.
   *
   * @param args The arguments to the XOR constraint.
   */
  public void assertXor(Lit... args) {
    validate(args);
    if (args.length == 0) {
      MonosatJNI.AssertXors(this.getSolverPtr(), getBuffer(0, 0), 0);
    } else if (args.length == 1) {
      assertTrue(args[0]); // If this the correct behaviour here?
    } else if (args.length == 2) {
      MonosatJNI.AssertXor(this.getSolverPtr(), args[0].l, args[1].l);
    } else {
      assertXor(Arrays.asList(args));
    }
  }

  /**
   * Assert that an odd number of args must hold. Trivially contradictory if args is empty.
   *
   * @param args The arguments to the XOR constraint.
   */
  public void assertXor(Collection<Lit> args) {
    validate(args);
    MonosatJNI.AssertXors(this.getSolverPtr(), this.getLitBuffer(args), args.size());
  }

  /**
   * Assert that an even number of args must hold.
   *
   * @param args The arguments to the XNOR constraint.
   */
  public void assertXnor(Lit... args) {
    validate(args);
    if (args.length == 0) {
      MonosatJNI.AssertXnors(this.getSolverPtr(), getBuffer(0, 0), 0);
    } else if (args.length == 1) {
      assertTrue(args[0]); // If this the correct behaviour here?
    } else if (args.length == 2) {
      MonosatJNI.AssertXnor(this.getSolverPtr(), args[0].l, args[1].l);
    } else {
      assertXnor(Arrays.asList(args));
    }
  }

  /**
   * Assert that an even number of args must hold.
   *
   * @param args The arguments to the XNOR constraint.
   */
  public void assertXnor(Collection<Lit> args) {
    validate(args);
    MonosatJNI.AssertXnors(this.getSolverPtr(), this.getLitBuffer(args), args.size());
  }

  /**
   * Assert that a and b must have the save value.
   *
   * @param a The first argument to make equal.
   * @param b The second argument to make equal.
   */
  public void assertEqual(Lit a, Lit b) {
    validate(a, b);
    assertXnor(a, b);
  }

  /**
   * Assert that a implies b: If a is true, then b must be true.
   *
   * @param a The pre-condition.
   * @param b The post-condition.
   */
  public void assertImplies(Lit a, Lit b) {
    validate(a, b);
    assertOr(a.not(), b);
  }

  /**
   * Assert that a implies that at least one element of args is true
   *
   * @param a The pre-condition.
   * @param args The post-condition disjunction.
   */
  public void assertImpliesOr(Lit a, Collection<Lit> args) {
    validate(a);
    validate(args);

      MonosatJNI.AssertImpliesOr(this.getSolverPtr(), a.toInt(), this.getLitBuffer(args), args.size());

  }
  /**
   * Assert that a implies that at least one element of args is true
   *
   * @param a The pre-condition.
   * @param args The post-condition disjunction.
   */
  public void assertImpliesOr(Lit a, Lit... args) {
    assertImpliesOr(a,Arrays.asList(args));
  }

  /**
   * Assert that a implies that all elements of args are true
   *
   * @param a The pre-condition.
   * @param args The post-condition disjunction.
   */
  public void assertImpliesAnd(Lit a, Collection<Lit> args) {
    validate(a);
    validate(args);
    if(args.size()==0){
      //do nothing
    }else {
      MonosatJNI.AssertImpliesAnd(this.getSolverPtr(), a.toInt(), this.getLitBuffer(args), args.size());
    }
  }

  /**
   * Assert that a implies that all elements of args are true
   *
   * @param a The pre-condition.
   * @param args The post-condition disjunction.
   */
  public void assertImpliesAnd(Lit a, Lit... args) {
    validate(a);
    validate(args);
    if(args.length==0){
      //do nothing
    }else{
      MonosatJNI.AssertImpliesAnd(this.getSolverPtr(), a.toInt(), this.getLitBuffer(args,0), args.length);
    }
  }

  /**
   * Create a new BitVector that will evaluate to the same value as 'then' if condition is true, and
   * to 'els' otherwise.
   *
   * @param condition The condition literal to test.
   * @param then The value of the returned BitVector if 'condition' is true.
   * @param els The value of the returned BitVector if 'condition' is false.
   * @return A new literal, equal to 'then' if condition is true, and equal to 'els' if condition is
   *     false.
   */
  public BitVector ite(Lit condition, BitVector then, BitVector els) {
    validate(condition);
    validate(then, els);
    assert (then.width() == els.width());
    BitVector result = new BitVector(this, then.width());
    MonosatJNI.bv_ite(this.getSolverPtr(), this.bvPtr, condition.l, then.id, els.id, result.id);
    return result;
  }

  // BitVector constructs

  /**
   * Return a new BitVector, equal to the bit-wise AND of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise AND of a and b.
   */
  public BitVector and(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_and(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Return a new BitVector, equal to the bit-wise OR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise OR of a and b.
   */
  public BitVector or(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_or(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Return a new BitVector, equal to the bit-wise NOT of 'a'.
   *
   * @param a The BitVector to bit-wise invert.
   * @return A bitvector of the same width as a, equal to the bit-wise NOT of a.
   */
  public BitVector not(BitVector a) {
    validate(a);
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_not(getSolverPtr(), bvPtr, a.id, result.id);
    return result;
  }

  /**
   * Return a new BitVector, equal to the bit-wise NAND of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise NAND of a and b.
   */
  public BitVector nand(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_nand(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Return a new BitVector, equal to the bit-wise NOR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise NOR of a and b.
   */
  public BitVector nor(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_nor(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Return a new BitVector, equal to the bit-wise XOR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise XOR of a and b.
   */
  public BitVector xor(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_xor(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Return a new BitVector, equal to the bit-wise XNOR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise XNOR of a and b.
   */
  public BitVector xnor(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_xnor(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Returns a BitVector that represents the non-wrapping two's complement addition of a and b. To
   * prevent wrapping, the solver will enforce that a+b<2^width.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A new BitVector with the same width as a, equal to a + b.
   */
  public BitVector add(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_addition(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Returns a BitVector that represents the non-wrapping two's complement subtraction of this and
   * other. To prevent wrapping, the solver will enforce that a-b>=0.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A new BitVector with the same width as a, equal to a - b.
   */
  public BitVector subtract(BitVector a, BitVector b) {
    validate(a, b);
    assert (a.width() == b.width());
    BitVector result = new BitVector(this, a.width());
    MonosatJNI.bv_subtraction(getSolverPtr(), bvPtr, a.id, b.id, result.id);
    return result;
  }

  /**
   * Assert that BitVector's a and b are equal to each other.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   */
  public void assertEqual(BitVector a, BitVector b) {
    validate(a, b);
    assertTrue(a.eq(b));
  }

  /**
   * Assert that BitVector a is equal to a constant.
   *
   * @param a The BitVector to be constrained.
   * @param constant Must be non-negative, and < 2^a.width()
   */
  public void assertEqual(BitVector a, long constant) {
    validate(a);
    assertTrue(a.eq(constant));
  }

  /**
   * Assert that BitVector a is equal to a constant.
   *
   * @param constant Must be non-negative, and < 2^a.width()
   * @param a The BitVector to be constrained.
   */
  public void assertEqual(long constant, BitVector a) {
    validate(a);
    assertTrue(a.eq(constant));
  }

  /**
   * Create a new BitVector, and assert that it is equal to the smallest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public BitVector min(Collection<BitVector> args) {
    validateBV(args);
    if (args.size() == 0) {
      throw new IllegalArgumentException("args must be non-empty");
    }
    int w = args.iterator().next().width();
    BitVector result = new BitVector(this, w);
    MonosatJNI.bv_min(
        this.getSolverPtr(), this.bvPtr, this.getBVBuffer(args, 0), args.size(), result.id);
    return result;
  }

  /**
   * Create a new BitVector, and assert that it is equal to the smallest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public BitVector min(BitVector... args) {
    return max(Arrays.asList(args));
  }

  /**
   * Create a new BitVector, and assert that it is equal to the largest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public BitVector max(Collection<BitVector> args) {
    validateBV(args);
    if (args.size() == 0) {
      throw new IllegalArgumentException("args must be non-empty");
    }
    int w = args.iterator().next().width();
    BitVector result = new BitVector(this, w);
    MonosatJNI.bv_min(
        this.getSolverPtr(), this.bvPtr, this.getBVBuffer(args, 0), args.size(), result.id);
    return result;
  }

  /**
   * Create a new BitVector, and assert that it is equal to the largest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public BitVector max(BitVector... args) {
    return max(Arrays.asList(args));
  }

  /**
   * Represents a value that is either true, false, or undefined. Only for internal use, for
   * interfacing with the native library.
   */
  protected enum LBool {
    // Don't change the order of these, as they must match the order of the l_bool enum defined in
    // Monosat.
    True,
    False,
    Undef;

    static {
      True.opt = Optional.of(true);
      False.opt = Optional.of(false);
      Undef.opt = Optional.empty();
    }

    /**
     * Optional value objects are created statically allocated for each LBool type, so that we don't
     * need to instantiate new Optional instances on the heap for toOpt() calls. This matters, as a
     * user may look up the values of millions of literals in common use cases.
     */
    @SuppressWarnings("OptionalUsedAsFieldOrParameterType")
    private Optional<Boolean> opt;

    /**
     * Convert an integer to an LBool.
     *
     * @param value The integer representation of the LBool.
     * @return The LBool constant corresponding to value.
     */
    public static LBool toLbool(int value) {
      return values()[value];
    }

    /**
     * Convert a Boolean to either LBool.True or LBool.False.
     *
     * @param value The Boolean to convert.
     * @return The LBool constant corresponding to value.
     */
    public static LBool fromBool(boolean value) {
      return values()[value ? 0 : 1];
    }

    /**
     * Convert this LBool into an Optional Boolean.
     *
     * @return An Optional representing the Boolean value of this LBool.
     */
    public Optional<Boolean> toOpt() {
      return opt;
    }
  }
}
