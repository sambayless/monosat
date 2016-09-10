package synoptic.algorithms;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.List;

import monosat.MonosatLibrary;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;

public class Monosat {
	public class FailedToSolveException extends Exception {
		public FailedToSolveException(String message) {
			super(message);
		}

		public FailedToSolveException() {
			super("Failed to solve instance (may be SAT or UNSAT)");
		}
	}

	private com.sun.jna.Pointer solver;
	private com.sun.jna.Pointer fsmTheory = null;
	private com.sun.jna.Pointer bvTheory;
	private ArrayList<Integer> opts = new ArrayList<Integer>();
	private ByteBuffer store = null;
	private ByteBuffer store2 = null;

	private ArrayList<com.sun.jna.Pointer> graphs = new ArrayList<com.sun.jna.Pointer>();

	public Monosat() {
		solver = MonosatLibrary.newSolver();
	}

	private void initFSMTheory() {
		if (fsmTheory == null) {
			fsmTheory = MonosatLibrary.initFSMTheory(solver);
		}
		assert (fsmTheory != null);
	}

	void disablePreprocessing() {
		MonosatLibrary.disablePreprocessing(solver);
	}

	public Monosat(String args) {
		// StringPointer ptr = new StringPointer(args);
		// solver = MonosatLibrary.newSolver_arg(ptr.byValue());
		System.out.println("Expected solver args " + args);
		Pointer ptr = new Memory(args.length() + 1);
		ptr.setString(0, args);
		solver = MonosatLibrary.newSolver_arg(ptr);
		bvTheory = MonosatLibrary.initBVTheory(solver);
	}

	public void setOutputFile(String filename) {
		// StringPointer ptr = new StringPointer(filename);
		System.out.println("Expected output file " + filename);
		Pointer m = new Memory(filename.length() + 1);

		m.setString(0, filename);
		MonosatLibrary.setOutputFile(solver, m);
	}

	// Convert a dimacs style, +/- non-zero signed literal into a Minisat style
	// unsigned lit
	private int intToLit(int l) {
		if (l == 0) {
			throw new RuntimeException("Literals must be non-zero");
		}
		assert (l != 0);
		assert (Math.abs(l) - 1 < nVars());
		int sign = l < 0 ? 1 : 0;
		return (Math.abs(l) - 1) * 2 + sign;
	}

	private int intToVar(int l) {
		if (l <= 0) {
			throw new RuntimeException(
					"Only positive literals can be converted into variables, but literal "
							+ l + " is not positivte");
		}
		assert (l != 0);
		assert (Math.abs(l) - 1 < nVars());
		return Math.abs(l) - 1;
	}

	private int litToInt(int l) {
		assert (l >= 0);
		int sign = (l & 1) > 0 ? -1 : 1;
		int r = ((l / 2) + 1) * sign;
		assert (intToLit(r) == l);
		return r;
	}

	private IntBuffer getIntBuf(int[] clause) {
		if (store == null || (store.asIntBuffer().capacity() < clause.length)) {
			// force store to at least double in size
			store = ByteBuffer.allocateDirect(clause.length * 4 * 2);
			store.order(ByteOrder.nativeOrder());
		}
		IntBuffer b = store.asIntBuffer();
		b.clear();
		b.put(clause);
		return b;
	}

	private IntBuffer getVarBuf(int[] clause) {
		if (store == null || (store.asIntBuffer().capacity() < clause.length)) {
			// force store to at least double in size
			store = ByteBuffer.allocateDirect(clause.length * 4 * 2);
			store.order(ByteOrder.nativeOrder());
		}
		IntBuffer b = store.asIntBuffer();
		b.clear();
		for (int i = 0; i < clause.length; i++) {
			b.put(intToVar(clause[i]));
		}
		return b;
	}

	private IntBuffer getLitBuf(int[] clause) {
		if (store == null || (store.asIntBuffer().capacity() < clause.length)) {
			// force store to at least double in size
			store = ByteBuffer.allocateDirect(clause.length * 4 * 2);
			store.order(ByteOrder.nativeOrder());
		}
		IntBuffer b = store.asIntBuffer();
		b.clear();
		for (int i = 0; i < clause.length; i++) {
			b.put(intToLit(clause[i]));
		}
		return b;
	}

	private IntBuffer getIntBuf(List<Integer> clause) {
		if (store == null || (store.asIntBuffer().capacity() < clause.size())) {
			// force store to at least double in size
			store = ByteBuffer.allocateDirect(clause.size() * 4 * 2);
			store.order(ByteOrder.nativeOrder());
		}
		IntBuffer b = store.asIntBuffer();
		b.clear();
		for (int l : clause) {
			b.put(l);
		}
		return b;
	}

	private IntBuffer getIntBuf2(List<Integer> clause) {
		if (store2 == null || (store2.asIntBuffer().capacity() < clause.size())) {
			// force store to at least double in size
			store2 = ByteBuffer.allocateDirect(clause.size() * 4 * 2);
			store2.order(ByteOrder.nativeOrder());
		}
		IntBuffer b = store2.asIntBuffer();
		b.clear();
		for (int l : clause) {
			b.put(l);
		}
		return b;
	}

	private IntBuffer getVarBuf(List<Integer> clause) {
		if (store == null || (store.asIntBuffer().capacity() < clause.size())) {
			// force store to at least double in size
			store = ByteBuffer.allocateDirect(clause.size() * 4 * 2);
			store.order(ByteOrder.nativeOrder());
		}
		IntBuffer b = store.asIntBuffer();
		b.clear();
		for (int l : clause) {
			b.put(intToVar(l));
		}
		return b;
	}

	private IntBuffer getLitBuf(List<Integer> clause) {
		if (store == null || (store.asIntBuffer().capacity() < clause.size())) {
			// force store to at least double in size
			store = ByteBuffer.allocateDirect(clause.size() * 4 * 2);
			store.order(ByteOrder.nativeOrder());
		}
		IntBuffer b = store.asIntBuffer();
		b.clear();
		for (int l : clause) {
			b.put(intToLit(l));
		}
		return b;
	}

	public boolean solve(int... assumptions) throws FailedToSolveException {

		int r;
		if (opts.size() == 0) {
			r = MonosatLibrary.solveAssumptionsLimited(solver,
					getLitBuf(assumptions), assumptions.length);
		} else {
			r = MonosatLibrary.solveAssumptionsLimited_MinBVs(solver,
					getLitBuf(assumptions), assumptions.length,
					getIntBuf2(opts), opts.size());
		}

		if (r == 0)
			return true;// instance is SAT
		else if (r == 1)
			return false;// instance is UNSAT (under assumptions)
		else
			throw new FailedToSolveException();// Could not solve instance under
												// resource limitations
	}

	public void backtrack() {
		MonosatLibrary.backtrack(solver);
	}

	public void addClause(List<Integer> lits) {
		MonosatLibrary.addClause(solver, getLitBuf(lits), lits.size());
	}

	public void addClause(int... lits) {
		// ByteBuffer b = getBuf(clause);
		// Pointer imagePointer = Native.getDirectBufferPointer(b);
		// IntByReference ref = new IntByReference(imagePointer);
		MonosatLibrary.addClause(solver, getLitBuf(lits), lits.length);
	}

	public void minimize_bv(int bv) {
		opts.add(bv);
	}

	public void maximize_bv(int bv) {
		int to_minimize = newBitvector_anon(bv_width(bv));
		bv_not(bv, to_minimize);
		opts.add(to_minimize);
	}

	public void clearOptimizationArgs() {
		opts.clear();
	}

	public void addClause(int a) {
		MonosatLibrary.addUnitClause(solver, intToLit(a));
	}

	public void addClause(int a, int b) {
		MonosatLibrary.addBinaryClause(solver, intToLit(a), intToLit(b));
	}

	public void addClause(int a, int b, int c) {
		MonosatLibrary.addTertiaryClause(solver, intToLit(a), intToLit(b),
				intToLit(c));
	}

	public int newLit() {
		return litToInt(MonosatLibrary.newVar(solver) * 2);
	}

	public int nVars() {
		return MonosatLibrary.nVars(solver);
	}

	public int newBitvector_const(int bvWidth, long constval) {
		assert (bvWidth > 0);
		assert (constval >= 0);// only unsigned values here
		return MonosatLibrary.newBitvector_const(solver, bvTheory, bvWidth,
				constval);
	}

	public int newBitvector_anon(int width) {
		return MonosatLibrary.newBitvector_anon(solver, bvTheory, width);
	}

	public int newBitvector(int width) {
		ArrayList<Integer> bits = new ArrayList<Integer>();
		for (int i = 0; i < width; i++) {
			bits.add(newLit());
		}
		return newBitvector(bits);
	}

	public int newBitvector(int[] bits) {
		return MonosatLibrary.newBitvector(solver, bvTheory, getVarBuf(bits),
				bits.length);
	}

	public int newBitvector(List<Integer> bits) {
		return MonosatLibrary.newBitvector(solver, bvTheory, getVarBuf(bits),
				bits.size());
	}

	public int newBitvector_comparison_const_lt(int bvID, long weight) {
		return litToInt(MonosatLibrary.newBVComparison_const_lt(solver,
				bvTheory, bvID, weight));
	}

	public int newBitvector_comparison_const_leq(int bvID, long weight) {

		return litToInt(MonosatLibrary.newBVComparison_const_leq(solver,
				bvTheory, bvID, weight));
	}

	public int newBitvector_comparison_const_gt(int bvID, long weight) {
		return litToInt(MonosatLibrary.newBVComparison_const_gt(solver,
				bvTheory, bvID, weight));
	}

	public int newBitvector_comparison_const_geq(int bvID, long weight) {
		return litToInt(MonosatLibrary.newBVComparison_const_geq(solver,
				bvTheory, bvID, weight));
	}

	public int newBitvector_comparison_bv_lt(int bvID, int bvIDto) {
		return litToInt(MonosatLibrary.newBVComparison_bv_lt(solver, bvTheory,
				bvID, bvIDto));
	}

	public int newBitvector_comparison_bv_leq(int bvID, int bvIDto) {
		return litToInt(MonosatLibrary.newBVComparison_bv_leq(solver, bvTheory,
				bvID, bvIDto));
	}

	public int newBitvector_comparison_bv_gt(int bvID, int bvIDto) {
		return litToInt(MonosatLibrary.newBVComparison_bv_gt(solver, bvTheory,
				bvID, bvIDto));
	}

	public int newBitvector_comparison_bv_geq(int bvID, int bvIDto) {
		return litToInt(MonosatLibrary.newBVComparison_bv_geq(solver, bvTheory,
				bvID, bvIDto));
	}

	public int bv_width(int bvID) {
		return MonosatLibrary.bv_width(solver, bvTheory, bvID);
	}

	public int bv_concat(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_concat(aID, bID, resultID);
		return resultID;
	}

	public void bv_concat(int aID, int bID, int resultID) {
		MonosatLibrary.bv_concat(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_not(int aID) {
		int resultID = newBitvector_anon(bv_width(aID));
		bv_not(aID, resultID);
		return resultID;
	}

	public void bv_not(int aID, int resultID) {
		MonosatLibrary.bv_not(solver, bvTheory, aID, resultID);
	}

	public int bv_and(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_and(aID, bID, resultID);
		return resultID;
	}

	public void bv_and(int aID, int bID, int resultID) {
		MonosatLibrary.bv_and(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_nand(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_nand(aID, bID, resultID);
		return resultID;
	}

	public void bv_nand(int aID, int bID, int resultID) {
		MonosatLibrary.bv_nand(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_or(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_or(aID, bID, resultID);
		return resultID;
	}

	public void bv_or(int aID, int bID, int resultID) {
		MonosatLibrary.bv_or(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_nor(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_nor(aID, bID, resultID);
		return resultID;
	}

	public void bv_nor(int aID, int bID, int resultID) {
		MonosatLibrary.bv_nor(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_xor(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_xor(aID, bID, resultID);
		return resultID;
	}

	public void bv_xor(int aID, int bID, int resultID) {
		MonosatLibrary.bv_xor(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_xnor(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_xnor(aID, bID, resultID);
		return resultID;
	}

	public void bv_xnor(int aID, int bID, int resultID) {
		MonosatLibrary.bv_xnor(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_ite(int condition_lit, int bvThenID, int bvElseID) {
		assert (bv_width(bvThenID) == bv_width(bvElseID));
		int resultID = newBitvector(bv_width(bvThenID));
		bv_ite(condition_lit, bvThenID, bvElseID, resultID);
		return resultID;
	}

	public void bv_ite(int conditionLit, int bvThenID, int bvElseID,
			int bvResultID) {
		MonosatLibrary.bv_ite(solver, bvTheory, conditionLit, bvThenID,
				bvElseID, bvResultID);
	}

	public int bv_addition(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_addition(aID, bID, resultID);
		return resultID;
	}

	public void bv_addition(int aID, int bID, int resultID) {
		MonosatLibrary.bv_addition(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_subtraction(int aID, int bID) {
		assert (bv_width(aID) == bv_width(bID));
		int resultID = newBitvector(bv_width(aID));
		bv_subtraction(aID, bID, resultID);
		return resultID;
	}

	public void bv_subtraction(int aID, int bID, int resultID) {
		MonosatLibrary.bv_subtraction(solver, bvTheory, aID, bID, resultID);
	}

	public int bv_slice(int aID, int lower, int upper) {
		assert (upper >= lower);
		int resultID = newBitvector(upper - lower);
		bv_slice(aID, lower, upper, resultID);
		return resultID;
	}

	public void bv_slice(int aID, int lower, int upper, int resultID) {
		MonosatLibrary.bv_slice(solver, bvTheory, aID, lower, upper, resultID);
	}

	public int bv_min(int... args) {
		assert (args.length > 0);
		int resultID = newBitvector(bv_width(args[0]));
		bv_min(resultID, args);
		return resultID;
	}

	public void bv_min(int resultID, int... args) {
		MonosatLibrary.bv_min(solver, bvTheory, getIntBuf(args), args.length,
				resultID);
	}

	public int bv_min(List<Integer> args) {
		assert (args.size() > 0);
		int resultID = newBitvector(bv_width(args.get(0)));
		bv_min(resultID, args);
		return resultID;
	}

	public void bv_min(int resultID, List<Integer> args) {
		MonosatLibrary.bv_min(solver, bvTheory, getIntBuf(args), args.size(),
				resultID);
	}

	public int bv_max(int... args) {
		assert (args.length > 0);
		int resultID = newBitvector(bv_width(args[0]));
		bv_max(resultID, args);
		return resultID;
	}

	public void bv_max(int resultID, int... args) {
		MonosatLibrary.bv_max(solver, bvTheory, getIntBuf(args), args.length,
				resultID);
	}

	public int bv_max(List<Integer> args) {
		assert (args.size() > 0);
		int resultID = newBitvector(bv_width(args.get(0)));
		bv_max(resultID, args);
		return resultID;
	}

	public void bv_max(int resultID, List<Integer> args) {
		MonosatLibrary.bv_max(solver, bvTheory, getIntBuf(args), args.size(),
				resultID);
	}

	public int bv_popcount(int... args) {
		assert (args.length > 0);
		int bvWidth = Integer.highestOneBit(args.length);
		int resultID = newBitvector_anon(bvWidth);
		bv_popcount(resultID, args);
		return resultID;
	}

	public void bv_popcount(int resultID, int... args) {
		MonosatLibrary.bv_popcount(solver, bvTheory, getLitBuf(args),
				args.length, resultID);
	}

	public int bv_popcount(List<Integer> args) {
		assert (args.size() > 0);
		int bvWidth = Integer.highestOneBit(args.size());
		int resultID = newBitvector_anon(bvWidth);
		bv_popcount(resultID, args);
		return resultID;
	}

	public void bv_popcount(int resultID, List<Integer> args) {
		MonosatLibrary.bv_popcount(solver, bvTheory, getLitBuf(args),
				args.size(), resultID);
	}

	public int newGraph() {
		Pointer g = MonosatLibrary.newGraph(solver);
		graphs.add(g);
		return graphs.size() - 1;
	}

	public int newNode(int graphID) {
		return MonosatLibrary.newNode(solver, graphs.get(graphID));
	}

	public int newEdge(int graphID, int from, int to) {
		return newEdge(graphID, from, to, 1);
	}

	public int newEdge(int graphID, int from, int to, long weight) {
		return litToInt(MonosatLibrary.newEdge(solver, graphs.get(graphID),
				from, to, weight));
	}

	public int reaches(int graphID, int from, int to) {
		return litToInt(MonosatLibrary.reaches(solver, graphs.get(graphID),
				from, to));
	}

	public int shortestPath_lt(int graphID, int from, int to, long lt) {
		return litToInt(MonosatLibrary.shortestPath_lt_const(solver,
				graphs.get(graphID), from, to, lt));
	}

	public int shortestPath_leq(int graphID, int from, int to, long leq) {
		return litToInt(MonosatLibrary.shortestPath_leq_const(solver,
				graphs.get(graphID), from, to, leq));
	}

	public int newFiniteStateMachine(int inputAlphabet) {
		return newFiniteStateMachine(inputAlphabet, 0);
	}

	public int newFiniteStateMachine(int inputAlphabet, int outputAlphabet) {
		initFSMTheory();
		return MonosatLibrary.newFSM(solver, fsmTheory, inputAlphabet,
				outputAlphabet);
	}

	public int newState(int fsmID) {
		initFSMTheory();
		return MonosatLibrary.newState(solver, fsmTheory, fsmID);
	}

	public int newTransition(int fsmID, int from, int to, int label) {
		return newTransition(fsmID, from, to, label, 0);
	}

	public int newTransition(int fsmID, int from, int to, int input, int output) {
		initFSMTheory();
		return litToInt(MonosatLibrary.newTransition(solver, fsmTheory, fsmID,
				from, to, input, output));
	}

	public int newString(List<Integer> string) {
		initFSMTheory();
		return MonosatLibrary.newString(solver, fsmTheory, getIntBuf(string),
				string.size());
	}

	public int newString(int... string) {
		initFSMTheory();
		return MonosatLibrary.newString(solver, fsmTheory, getIntBuf(string),
				string.length);
	}

	public int acceptsString(int fsmID, int startState, int acceptState,
			int stringID) {
		initFSMTheory();
		return litToInt(MonosatLibrary.fsmAcceptsString(solver, fsmTheory,
				fsmID, startState, acceptState, stringID));
	}

	public void atMostOne(List<Integer> lits) {
		ArrayList<Integer> vars = new ArrayList<Integer>();
		for (int l : lits) {
			assert (l != 0);
			int v;
			if (l < 0) {
				v = intToVar(newLit());
			} else {
				v = intToVar(l);
			}
			vars.add(v);
		}
		assert (vars.size() == lits.size());
		MonosatLibrary.at_most_one(solver, getIntBuf(vars), vars.size());
	}

	public void atMostOne(int... lits) {
		ArrayList<Integer> vars = new ArrayList<Integer>();
		for (int l : lits) {
			assert (l != 0);
			int v;
			if (l < 0) {
				v = intToVar(newLit());
			} else {
				v = intToVar(l);
			}
			vars.add(v);
		}
		assert (vars.size() == lits.length);
		MonosatLibrary.at_most_one(solver, getIntBuf(vars), vars.size());
	}

	// Get the assignment of a variable.
	// Must be called after solve(), and before any further constraints are
	// added to the solver.
	// defaults to false if the literal is unassigned
	public boolean value(int lit) {
		int v = MonosatLibrary.getModel_Literal(solver, intToLit(lit));
		if (v == 0)
			return true;
		else
			return false;
	}

	public boolean isAssigned(int lit) {
		int v = MonosatLibrary.getModel_Literal(solver, intToLit(lit));
		if (v == 2)
			return false;
		else
			return true;
	}

	// Get an assignment to a bitvector in the model. The model may find a range
	// of satisfying assignments to the bitvector; this function returns the
	// minimum satisfying value found the solver
	public long bv_value(int bvID) {
		return MonosatLibrary.getModel_BV(solver, bvTheory, bvID, false);
	}

	// Get an assignment to a bitvector in the model. The model may find a range
	// of satisfying assignments to the bitvector;
	// If getMaximumValue is true, this function returns the maximum satisfying
	// assignment to the bitvector in the model; else it returns the smallest.
	public long bv_value(int bvID, boolean getMaximumValue) {
		return MonosatLibrary.getModel_BV(solver, bvTheory, bvID,
				getMaximumValue);
	}
}