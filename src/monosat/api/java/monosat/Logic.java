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

import java.util.*;
import java.util.logging.Level;

import static java.util.Collections.emptyIterator;

/**
 * Logic provides static accessors for common logic operations. These static methods form a
 * light-weight domain specific language.
 *
 * <p>The expected usage is for all methods of Logic to be statically imported:
 *
 * <blockquote>
 *
 * <pre>{@code
 * import static monosat.Logic.*
 * }</pre>
 *
 * </blockquote>
 *
 * Statically importing the methods of this class makes standard Boolean logic operators such as
 * and(),or(), not(), as well as constants True and False, available, allowing for concise formula
 * descriptions.
 */
public final class Logic {
  /** True literal, shared among all solver instances. */
  public static final Lit True =
      Lit.True; // Is there any possibility that this won't be constructed in the proper order?
  /** False literal, shared among all solver instances. */
  public static final Lit False = Lit.False;
  // make these constants available to users who import monosat.Logic.*
  public static final Comparison EQ = Comparison.EQ;
  public static final Comparison NEQ = Comparison.NEQ;
  public static final Comparison GEQ = Comparison.GEQ;
  public static final Comparison GT = Comparison.GT;
  public static final Comparison LEQ = Comparison.LEQ;
  public static final Comparison LT = Comparison.LT;
  /** If true, log a warning the first time a global contradiction occurs. */
  private static boolean warn_contradictions = true;
  /**
   * If true, throw an exception when global trivial contradictions occur, instead of logging a
   * warning.
   */
  private static boolean throw_contradictions = false;
  /** Private constructor, to prevent instances of logic from being constructed */
  private Logic() {}

  // Don't import Lit.Undef or Lit.Error

  /**
   * Helper method to retrieve the solver from an array of literals.
   *
   * @param args An array of literals.
   * @return The solver of one of the literals, if at least one of the literals has a solver.
   */
  private static Solver getSolver(Lit... args) {
    for (Lit l : args) {
      if (l == null || l == Lit.Error || l == Lit.Undef) {
        throw new IllegalArgumentException("Invalid literal: " + String.valueOf(l));
      }
      if (l.solver != null) {
        return l.solver;
      }
    }
    return null;
  }



  /**
   * Helper method to retrieve the solver from an array of literals.
   *
   * @param args A collection of literals.
   * @return The solver of one of the literals, if at least one of the literals has a solver.
   */
  private static Solver getSolver(Iterable<Lit> args) {
    for (Lit l : args) {
      if (l == null || l == Lit.Error || l == Lit.Undef) {
        throw new IllegalArgumentException("Invalid literal: " + String.valueOf(l));
      }
      if (l.solver != null) {
        return l.solver;
      }
    }
    return null;
  }
  /**
   * Helper method to retrieve the solver from an array of literals.
   *
   * @param args A collection of literals.
   * @return The solver of one of the literals, if at least one of the literals has a solver.
   */
  private static Solver getSolver(Iterable<Lit> args, Lit... lits) {
    Solver s = getSolver(args);
    if (s== null){
      s = getSolver(lits);
    }
    return s;
  }
  /**
   * Helper method to retrieve the solver from an array of literals.
   *
   * @param args A collection of literals.
   * @return The solver of one of the literals, if at least one of the literals has a solver.
   */
  private static Solver getSolver(Lit[] args, Lit... lits) {
    Solver s = getSolver(args);
    if (s== null){
      s = getSolver(lits);
    }
    return s;
  }
  /**
   * By default, certain trivial contradictions (such as AssertFalse(True), or AssertOr() with empty
   * arguments) will be caught, and log a warning. This is useful, as such contradictions are almost
   * always unintentional. Use this method to disable warnings on trivial contradictions.
   */
  public static synchronized void allowContradictions() {
    warn_contradictions = false;
    throw_contradictions = false;
  }

  /**
   * By default, certain trivial contradictions (such as AssertFalse(True), or AssertOr() with empty
   * arguments) will be caught, and log a warning. Use this method to instead throw exceptions on
   * trivial contradictions.
   */
  public static synchronized void disallowContradictions() {
    throw_contradictions = true;
  }

  private static synchronized void contradiction(Lit a) {
    // this global contradiction on True/False makes all existing solvers UNSAT
    for (Solver s : Solver.solvers.keySet()) {
      s.assertFalse(Lit.True);
    }
    if (throw_contradictions) {
      if (a == Lit.False) {
        throw new TrivialContradictionException(
            "Constant False asserted to be True (which is almost always an error).");
      } else if (a == Lit.True) {
        throw new TrivialContradictionException(
            "Constant True asserted to be False (which is almost always an error).");
      }
    }
    if (warn_contradictions) {
      warn_contradictions = false;
      if (a == Lit.False) {
        Solver.log.log(
            Level.WARNING, "Constant False asserted to be True (which is almost always an error).");
      } else if (a == Lit.True) {
        Solver.log.log(
            Level.WARNING, "Constant True asserted to be False (which is almost always an error).");
      }
    }
  }

  /**
   * Helper method, called when a trivial global contradiction is detected. Makes all existing SAT
   * solvers UNSAT, and possibly logs a warning or throws an exception.
   */
  private static synchronized void contradiction() {
    // this global contradiction on True/False makes all existing solvers UNSAT
    for (Solver s : Solver.solvers.keySet()) {
      s.assertFalse(Lit.True);
    }
    if (throw_contradictions) {
      throw new TrivialContradictionException(
          "Statically UNSAT assertion (which is almost always an error).");
    }
    if (warn_contradictions) {
      warn_contradictions = false;
      Solver.log.log(
          Level.WARNING, "Statically UNSAT assertion (which is almost always an error).");
    }
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
  public static Lit not(Lit a) {
    return a.not();
  }

  /**
   * Create a new literal that evaluates to true if a is false, or if a is true and b is true.
   *
   * @param a The pre-condition.
   * @param b The post-condition.
   * @return A new literal representing an implication gate.
   */
  public static Lit implies(Lit a, Lit b) {
    Solver solver = getSolver(a, b);
    if (solver != null) {
      return solver.implies(a, b);
    } else {
      if (a == Lit.True && b == Lit.False) {
        return Lit.False;
      } else {
        return Lit.True;
      }
    }
  }
  /**
   *  Create a new literal that evaluates to true if a is false, or if a is true and at least one element of args is
   *  true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public static Lit impliesOr(Lit a, Collection<Lit> args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      return solver.impliesOr(a, args);
    } else {
      if (a == Lit.True) {
        boolean allFalse = true;
        for(Lit l:args){
          if(!l.isConstFalse()){
            allFalse = false;
          }
        }
        if(!allFalse || args.size()==0) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
  }

  /**
   *  Create a new literal that evaluates to true if a is false, or if a is true and at least one element of args is
   *  true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public static Lit impliesOr(Lit a, Lit... args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      return solver.impliesOr(a, args);
    } else {
      if (a == Lit.True) {
        boolean allFalse = true;
        for(Lit l:args){
          if(!l.isConstFalse()){
            allFalse = false;
          }
        }
        if(!allFalse || args.length==0) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
  }

  /**
   *  Create a new literal that evaluates to true if a is false, or if a is true and all literals in args
   *  are true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public static Lit impliesAnd(Lit a, Collection<Lit> args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      return solver.impliesAnd(a, args);
    } else {
      if (a == Lit.True) {
        boolean anyFalse = false;
        for(Lit l:args){
          if(l.isConstFalse()){
            anyFalse = true;
          }
        }
        if(anyFalse) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
  }
  /**
   *  Create a new literal that evaluates to true if a is false, or if a is true and all literals in args
   *  are true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   * @return A new literal representing an implication gate.
   */
  public static Lit impliesAnd(Lit a, Lit... args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      return solver.impliesAnd(a, args);
    } else {
      if (a == Lit.True) {
        boolean anyFalse = false;
        for(Lit l:args){
          if(l.isConstFalse()){
            anyFalse = true;
          }
        }
        if(anyFalse) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
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
  public static Lit ite(Lit condition, Lit then, Lit els) {
    Solver solver = getSolver(condition, then, els);
    if (solver != null) {
      return solver.ite(condition, then, els);
    } else {
      if (condition == Lit.True) {
        return then;
      } else {
        return els;
      }
    }
  }

  /**
   * Create a new literal that is true if all the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the And gate.
   * @return A new literal, true if all the literals in args are true.
   */
  public static Lit and(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.and(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
  }

  /**
   * Create a new literal that is true if all the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the And gate.
   * @return A new literal, true if all the literals in args are true.
   */
  public static Lit and(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.and(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
  }
  /**
   * Create a new literal that is true if any of the literals in args are true. If args is empty,
   * return Lit.False.
   *
   * @param args The arguments to the Or gate.
   * @return A new literal, true if any of the literals in args are true.
   */
  public static Lit or(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.or(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.True) {
          return Lit.True;
        }
      }
      return Lit.False;
    }
  }
  /**
   * Create a new literal that is true if any of the literals in args are true. If args is empty,
   * return Lit.False.
   *
   * @param args The arguments to the Or gate.
   * @return A new literal, true if any of the literals in args are true.
   */
  public static Lit or(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.or(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.True) {
          return Lit.True;
        }
      }
      return Lit.False;
    }
  }

  /**
   * Create a new literal that is true if not all of the literals in args are true. If args is
   * empty, return Lit.False.
   *
   * @param args The arguments to the Nand gate.
   * @return A new literal, true if not all the literals in args are true.
   */
  public static Lit nand(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.nand(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          return Lit.True;
        }
      }
      return Lit.False;
    }
  }

  /**
   * Create a new literal that is true if not all of the literals in args are true. If args is
   * empty, return Lit.False.
   *
   * @param args The arguments to the Nand gate.
   * @return A new literal, true if not all the literals in args are true.
   */
  public static Lit nand(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.nand(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          return Lit.True;
        }
      }
      return Lit.False;
    }
  }

  /**
   * Create a new literal that is true if none of the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the Nor gate.
   * @return A new literal, true if none of the literals in args are true.
   */
  public static Lit nor(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.nor(args);
    } else {
      for (Lit l : args) {
        if (l != Lit.False) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
  }

  /**
   * Create a new literal that is true if none of the literals in args are true. If args is empty,
   * return Lit.True.
   *
   * @param args The arguments to the Nor gate.
   * @return A new literal, true if none of the literals in args are true.
   */
  public static Lit nor(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.nor(args);
    } else {
      for (Lit l : args) {
        if (l != Lit.False) {
          return Lit.False;
        }
      }
      return Lit.True;
    }
  }

  /**
   * Create a new literal that is true if an odd number of the literals in args are true. If args is
   * empty, return Lit.True.
   *
   * @param args The arguments to the Xor gate.
   * @return A new literal, true if an odd number of the literals in args are true.
   */
  public static Lit xor(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.xor(args);
    } else {
      // XOR is defined to be true if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 1) {
        return Lit.True;
      } else {
        return Lit.False;
      }
    }
  }

  /**
   * Create a new literal that is true if an odd number of the literals in args are true. If args is
   * empty, return Lit.True.
   *
   * @param args The arguments to the Xor gate.
   * @return A new literal, true if an odd number of the literals in args are true.
   */
  public static Lit xor(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.xor(args);
    } else {
      // XOR is defined to be true if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 1) {
        return Lit.True;
      } else {
        return Lit.False;
      }
    }
  }

  /**
   * Create a new literal that is true if an even number of the literals in args are true. If args
   * is empty, return Lit.True.
   *
   * @param args The arguments to the Xnor gate.
   * @return A new literal, true if an even number of the literals in args are true.
   */
  public static Lit xnor(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.xnor(args);
    } else {
      // XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 0) {
        return Lit.True;
      } else {
        return Lit.False;
      }
    }
  }

  /**
   * Create a new literal that is true if an even number of the literals in args are true. If args
   * is empty, return Lit.True.
   *
   * @param args The arguments to the Xnor gate.
   * @return A new literal, true if an even number of the literals in args are true.
   */
  public static Lit xnor(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      return solver.xnor(args);
    } else {
      // XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 0) {
        return Lit.True;
      } else {
        return Lit.False;
      }
    }
  }

  /**
   * Create a new literal, that is true if a and b both evaluate to the same truth value.
   *
   * @param a The first literal to test equality.
   * @param b The second literal to test equality.
   * @return A literal that is true if a and b both evaluate to the same truth value.
   */
  public static Lit equal(Lit a, Lit b) {
    return xnor(a, b);
  }
  // Assertion forms.

  /**
   * Assert that Lit a must be true in the solver.
   *
   * @param a The literal to assert true.
   */
  public static void assertTrue(Lit a) {
    Solver solver = getSolver(a);
    if (solver != null) {
      solver.assertTrue(a);
    } else {
      if (a == Lit.False) {
        contradiction(a);
      } else {
        // do nothing
      }
    }
  }

  /**
   * Assert that Lit a must be false in the solver.
   *
   * @param a The literal to assert false.
   */
  public static void assertFalse(Lit a) {
    Solver solver = getSolver(a);
    if (solver != null) {
      solver.assertFalse(a);
    } else {
      if (a == Lit.True) {
        contradiction(a);
      } else {
        // do nothing
      }
    }
  }

  /**
   * Assert that all of the literals in args are true. Trivially satisfied if args is empty.
   * Equivalent to assertAnd(args).
   *
   * @param args The arguments to assert true.
   */
  public static void assertTrue(Lit... args) {
    assertAnd(args);
  }

  /**
   * Assert that all of the literals in args are true. Trivially satisfied if args is empty.
   * Equivalent to assertAnd(args).
   *
   * @param args The arguments to assert true.
   */
  public static void assertTrue(Collection<Lit> args) {
    assertAnd(args);
  }

  /**
   * Assert that all of the literals in args are false. Trivially satisfied if args is empty.
   * Equivalent to assertNor(args).
   *
   * @param args The arguments to assert true.
   */
  public static void assertFalse(Lit... args) {
    assertNor(args);
  }

  /**
   * Assert that all of the literals in args are false. Trivially satisfied if args is empty.
   * Equivalent to assertNor(args).
   *
   * @param args The arguments to assert true.
   */
  public static void assertFalse(Collection<Lit> args) {
    assertNor(args);
  }

  /**
   * Assert that all of the literals in args are true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert true.
   */
  public static void assertAnd(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertAnd(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          contradiction(l);
        }
      }
      // do nothing
    }
  }
  /**
   * Assert that all of the literals in args are true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert true.
   */
  public static void assertAnd(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertAnd(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          contradiction(l);
        }
      }
      // do nothing
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
  public static void assertOr(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertOr(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.True) {
          return;
        }
      }
      contradiction();
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
  public static void assertOr(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertOr(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.True) {
          return;
        }
      }
      contradiction();
    }
  }

  /**
   * Assert that at least one of the given args must be false. Trivially contradictory if args is
   * empty.
   *
   * @param args The arguments to assert at least one must be false.
   */
  public static void assertNand(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertNand(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          return;
        }
      }
      contradiction();
    }
  }

  /**
   * Assert that at least one of the given args must be false. Trivially contradictory if args is
   * empty.
   *
   * @param args The arguments to assert at least one must be false.
   */
  public static void assertNand(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertNand(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.False) {
          return;
        }
      }
      contradiction();
    }
  }

  /**
   * Assert that none of the given args may be true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert false.
   */
  public static void assertNor(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertNor(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.True) {
          contradiction(l);
        }
      }
    }
  }

  /**
   * Assert that none of the given args may be true. Trivially satisfied if args is empty.
   *
   * @param args The arguments to assert false.
   */
  public static void assertNor(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertNor(args);
    } else {
      for (Lit l : args) {
        if (l == Lit.True) {
          contradiction(l);
        }
      }
    }
  }

  /**
   * Assert that an odd number of args must hold. Trivially contradictory if args is empty.
   *
   * @param args The arguments to the XOR constraint.
   */
  public static void assertXor(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertXor(args);
    } else {
      // XOR is defined to be true if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 1) {
        // do nothing
      } else {
        contradiction();
      }
    }
  }

  /**
   * Assert that an odd number of args must hold. Trivially contradictory if args is empty.
   *
   * @param args The arguments to the XOR constraint.
   */
  public static void assertXor(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertXor(args);
    } else {
      // XOR is defined to be true if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 1) {
        // do nothing
      } else {
        contradiction();
      }
    }
  }

  /**
   * Assert that an even number of args must hold.
   *
   * @param args The arguments to the XNOR constraint.
   */
  public static void assertXnor(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertXnor(args);
    } else {
      // XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 1) {
        contradiction();
      } else {
        // do nothing
      }
    }
  }

  /**
   * Assert that an even number of args must hold.
   *
   * @param args The arguments to the XNOR constraint.
   */
  public static void assertXnor(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertXnor(args);
    } else {
      // XNOR is defined to be false if an odd number of its arguments are true, and false otherwise
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount % 2 == 1) {
        contradiction();
      } else {
        // do nothing
      }
    }
  }

  /**
   * Assert that a and b must have the save value.
   *
   * @param a The first argument to make equal.
   * @param b The second argument to make equal.
   */
  public static void assertEqual(Lit a, Lit b) {
    Solver solver = getSolver(a, b);
    if (solver != null) {
      solver.assertEqual(a, b);
    } else {
      if (a != b) {
        contradiction();
      }
    }
  }

  /**
   * Assert that a implies b: If a is true, then b must be true.
   *
   * @param a The pre-condition.
   * @param b The post-condition.
   */
  public static void assertImplies(Lit a, Lit b) {
    Solver solver = getSolver(a, b);
    if (solver != null) {
      solver.assertImplies(a, b);
    } else {
      if (a == Lit.True && b == Lit.False) {
        contradiction();
      }
    }
  }

  /**
   * Assert that a implies the disjunction of args: If a is true, then at least one element of args must be true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   */
  public static void assertImpliesOr(Lit a, Collection<Lit> args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      solver.assertImpliesOr(a, args);
    } else {
      if (a == Lit.True) {
        boolean allFalse = true;
        for(Lit l:args){
          if(!l.isConstFalse()){
            allFalse = false;
          }
        }
        if(allFalse || args.size()==0) {
          contradiction();
        }
      }
    }
  }
  /**
   * Assert that a implies the disjunction of args: If a is true, then at least one element of args must be true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   */
  public static void assertImpliesOr(Lit a, Lit... args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      solver.assertImpliesOr(a, args);
    } else {
      if (a == Lit.True) {
        boolean allFalse = true;
        for(Lit l:args){
          if(!l.isConstFalse()){
            allFalse = false;
          }
        }
        if(allFalse || args.length==0) {
          contradiction();
        }
      }
    }
  }

  /**
   * Assert that a implies the conjunction of args: If a is true, then all elements elements of args must be true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   */
  public static void assertImpliesAnd(Lit a, Collection<Lit> args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      solver.assertImpliesAnd(a, args);
    } else {
      if (a == Lit.True) {
        boolean anyFalse = false;
        for(Lit l:args){
          if(l.isConstFalse()){
            anyFalse = true;
          }
        }
        if(anyFalse) {
          contradiction();
        }
      }
    }
  }
  /**
   * Assert that a implies the conjunction of args: If a is true, then all elements elements of args must be true.
   *
   * @param a The pre-condition.
   * @param args The post-condition.
   */
  public static void assertImpliesAnd(Lit a, Lit... args) {
    Solver solver = getSolver(args,a);
    if (solver != null) {
      solver.assertImpliesAnd(a, args);
    } else {
      if (a == Lit.True) {
        boolean anyFalse = false;
        for(Lit l:args){
          if(l.isConstFalse()){
            anyFalse = true;
          }
        }
        if(anyFalse) {
          contradiction();
        }
      }
    }
  }

  // Bitvector constructs
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
  public static BitVector ite(Lit condition, BitVector then, BitVector els) {
    return then.getSolver().ite(condition, then, els);
  }

  /**
   * Return a new BitVector, equal to the bit-wise AND of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise AND of a and b.
   */
  public static BitVector and(BitVector a, BitVector b) {
    return a.getSolver().and(a, b);
  }
  /**
   * Return a new BitVector, equal to the bit-wise OR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise OR of a and b.
   */
  public static BitVector or(BitVector a, BitVector b) {
    return a.getSolver().or(a, b);
  }
  /**
   * Return a new BitVector, equal to the bit-wise NOT of 'a'.
   *
   * @param a The BitVector to bit-wise invert.
   * @return A bitvector of the same width as a, equal to the bit-wise NOT of a.
   */
  public static BitVector not(BitVector a) {
    return a.getSolver().not(a);
  }
  /**
   * Return a new BitVector, equal to the bit-wise NAND of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise NAND of a and b.
   */
  public static BitVector nand(BitVector a, BitVector b) {
    return a.getSolver().nand(a, b);
  }
  /**
   * Return a new BitVector, equal to the bit-wise NOR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise NOR of a and b.
   */
  public static BitVector nor(BitVector a, BitVector b) {
    return a.getSolver().nor(a, b);
  }

  /**
   * Return a new BitVector, equal to the bit-wise XOR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise XOR of a and b.
   */
  public static BitVector xor(BitVector a, BitVector b) {
    return a.getSolver().xor(a, b);
  }

  /**
   * Return a new BitVector, equal to the bit-wise XNOR of 'a' and 'b'.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A bitvector of the same width as a, equal to the bit-wise XNOR of a and b.
   */
  public static BitVector xnor(BitVector a, BitVector b) {
    return a.getSolver().xnor(a, b);
  }

  /**
   * Returns a Bitvector that represents the non-wrapping two's complement addition of a and b. To
   * prevent wrapping, the solver will enforce that a+b<2^width.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A new Bitvector with the same width as a, equal to a + b.
   */
  public static BitVector add(BitVector a, BitVector b) {
    return a.getSolver().add(a, b);
  }
  /**
   * Returns a Bitvector that represents the non-wrapping two's complement subtraction of this and
   * other. To prevent wrapping, the solver will enforce that a-b>=0.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   * @return A new Bitvector with the same width as a, equal to a - b.
   */
  public static BitVector subtract(BitVector a, BitVector b) {
    return a.getSolver().subtract(a, b);
  }
  /**
   * Assert that BitVector's a and b are equal to each other.
   *
   * @param a The first argument. Must have the same bitwidth as 'b'.
   * @param b The second argument. Must have the same bitwidth as 'a'.
   */
  public static void assertEqual(BitVector a, BitVector b) {
    a.getSolver().assertEqual(a, b);
  }

  /**
   * Assert that BitVector's a is equal to a constant.
   *
   * @param a The BitVector to be constrained.
   * @param constant Must be non-negative, and < 2^a.width()
   */
  public static void assertEqual(BitVector a, long constant) {
    a.getSolver().assertEqual(a, constant);
  }
  /**
   * Assert that BitVector's a is equal to a constant.
   *
   * @param constant Must be non-negative, and < 2^a.width()
   * @param a The BitVector to be constrained.
   */
  public static void assertEqual(long constant, BitVector a) {
    a.getSolver().assertEqual(constant, a);
  }
  /**
   * Create a new BitVector, and assert that it is equal to the smallest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public static BitVector min(Collection<BitVector> args) {
    if (args.size() > 0) {
      return args.iterator().next().getSolver().min(args);
    } else {
      throw new IllegalArgumentException("min requires at least one argument");
    }
  }
  /**
   * Create a new BitVector, and assert that it is equal to the smallest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public static BitVector min(BitVector... args) {
    return min(Arrays.asList(args));
  }
  /**
   * Create a new BitVector, and assert that it is equal to the largest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public static BitVector max(Collection<BitVector> args) {
    if (args.size() > 0) {
      return args.iterator().next().getSolver().max(args);
    } else {
      throw new IllegalArgumentException("max requires at least one argument");
    }
  }
  /**
   * Create a new BitVector, and assert that it is equal to the largest element of args. Each
   * argument in args must have the same bitwidth.
   *
   * @param args A non-empty array of BitVectors. Must all have the same bitwidth.
   * @return A new BitVector, constrained to be equal to the smallest BitVector from args. Will have
   *     the same bitwidth as the elements of args.
   */
  public static BitVector max(BitVector... args) {
    return max(Arrays.asList(args));
  }

  /**
   * Enforce that at most one of the specified literals may be true. If 1 or fewer arguments are
   * given, has no effect. If exactly 2 arguments are given, this is the same as:
   * assertOr(args[0],args[1])
   *
   * @param args An array of Lits, at most one of which can be true.
   */
  public static void assertAtMostOne(Lit... args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertAtMostOne(args);
    } else {
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount > 1) contradiction();
    }
  }

  /**
   * Enforce that at most one of the specified literals may be true. If 1 or fewer arguments are
   * given, has no effect. If exactly 2 arguments are given, this is the same as:
   * assertOr(args[0],args[1])
   *
   * @param args A collection of Lits, at most one of which can be true.
   */
  public static void assertAtMostOne(Collection<Lit> args) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertAtMostOne(args);
    } else {
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      if (trueCount > 1) contradiction();
    }
  }

  /**
   * Enforce a pseudo-Boolean constraint. The number of true literals from among args must satisfy
   * comparison c relative to compare to.
   *
   * <p>For example, if c is Comparison.LEQ, and compareTo is 3, then at most 3 literals from args
   * may be true.
   *
   * @param args A collection of args, whose sum will be compared.
   * @param c The comparison operation to perform.
   * @param compareTo The constant to compare the sum of the true lits in args to.
   */
  public static void assertPB(Collection<Lit> args, Comparison c, int compareTo) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertPB(args, c, compareTo);
    } else {
      int trueCount = 0;
      for (Lit l : args) {
        if (l == Lit.True) {
          trueCount++;
        }
      }
      // statically check if the comparison holds
      if (!c.compare(trueCount, compareTo)) {
        contradiction();
      }
    }
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
  public static void assertPB(Collection<Lit> args, Collection<Integer> weights, Comparison c, int compareTo) {
    Solver solver = getSolver(args);
    if (solver != null) {
      solver.assertPB(args, c, compareTo);
    } else {
      int trueCount = 0;
      Iterator<Lit> it = args.iterator();
      Iterator<Integer> wit = weights!=null ? weights.iterator():Collections.emptyIterator();
      while (it.hasNext()) {
        Lit l = it.next();
        int weight = 1;
        if (wit.hasNext()) {
          weight = wit.next();
        }
        if (l == Lit.True) {
          trueCount += weight;
        }
      }
      // statically check if the comparison holds
      if (!c.compare(trueCount, compareTo)) {
        contradiction();
      }
    }
  }
}
