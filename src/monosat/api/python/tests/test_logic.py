import unittest
import warnings

import monosat


class TestLogic(unittest.TestCase):
    def test_newLit(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        self.assertNotEqual(a, monosat.logic.true())

        b = monosat.logic.Var()
        self.assertNotEqual(a, b)
        self.assertNotEqual(a.lit, b.lit)

        self.assertTrue(monosat.solver.Solve(a))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(b)))
        self.assertTrue(monosat.solver.Solve(a))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), b))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), a))
        self.assertTrue(monosat.solver.Solve(a))
        self.assertTrue(monosat.solver.Solve(b))

    def test_solve(self):
        monosat.Monosat().newSolver()
        self.assertTrue(monosat.solver.Solve())
        self.assertTrue(monosat.solver.Solve())
        monosat.logic.AssertTrue(monosat.logic.true())
        self.assertTrue(monosat.solver.Solve())
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            monosat.logic.AssertFalse(monosat.logic.true())
            self.assertTrue("Asserted a trivial contradiction" in str(w[-1].message))

        self.assertFalse(monosat.solver.Solve())
        monosat.Monosat().newSolver()
        self.assertTrue(monosat.solver.Solve())

    def test_deprecated_init(self):
        monosat.Monosat().newSolver()
        # tests the deprecated init function
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            monosat.Monosat().init("-decide-theories")
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertTrue("deprecated" in str(w[-1].message))

        self.assertTrue(monosat.solver.Solve())
        self.assertTrue(monosat.solver.Solve())
        monosat.logic.AssertTrue(monosat.logic.true())
        self.assertTrue(monosat.solver.Solve())
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            monosat.logic.AssertFalse(monosat.logic.true())
            self.assertTrue("Asserted a trivial contradiction" in str(w[-1].message))

        self.assertFalse(monosat.solver.Solve())
        monosat.Monosat().newSolver()
        self.assertTrue(monosat.solver.Solve())

    def test_testSolveAssumps(self):
        monosat.Monosat().newSolver()
        self.assertTrue(monosat.solver.Solve())
        self.assertTrue(monosat.solver.Solve())
        self.assertTrue(monosat.solver.Solve(monosat.logic.true()))
        self.assertFalse(monosat.solver.Solve(monosat.logic.false()))
        self.assertFalse(monosat.solver.Solve(monosat.logic.true(), monosat.logic.false()))
        self.assertTrue(monosat.solver.Solve())

    def test_ite(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        result = monosat.logic.Ite(c, a, b)
        self.assertTrue(monosat.solver.Solve(c, a, monosat.logic.Not(b), result))
        self.assertFalse(monosat.solver.Solve(c, a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(c, monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(result)))

        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(c), a, monosat.logic.Not(b), result))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(c), a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(c), monosat.logic.Not(a), b, result))

    def test_not(self):
        self.assertTrue(monosat.logic.Not(monosat.logic.true()), monosat.logic.false())
        a = monosat.logic.Var()
        self.assertNotEqual(monosat.logic.Not(a), a)
        self.assertTrue(monosat.logic.Not(a), a.Not())

    def test_and(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        result = monosat.logic.And(a, b)
        self.assertTrue(monosat.solver.Solve(a, b, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(a, b, result))

    def test_ands(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        result = monosat.logic.And(a, b, c)
        self.assertTrue(monosat.solver.Solve(a, b, c, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), c, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(a, b, c, result))

    def test_or(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        result = monosat.logic.Or(a, b)
        self.assertTrue(monosat.solver.Solve(a, b, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), result))
        self.assertTrue(monosat.solver.Solve(a, b, result))

    def test_ors(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        result = monosat.logic.Or(a, b, c)
        self.assertTrue(monosat.solver.Solve(a, b, c, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), c, result))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(c), result))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c, result))
        self.assertTrue(monosat.solver.Solve(a, b, c, result))

    def test_nand(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        result = monosat.logic.Nand(a, b)
        self.assertFalse(monosat.solver.Solve(a, b, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(a, b, result))

    def test_nands(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        result = monosat.logic.Nand(a, b, c)
        self.assertFalse(monosat.solver.Solve(a, b, c, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), c, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(a, b, c, result))

    def test_nor(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        result = monosat.logic.Nor(a, b)
        self.assertFalse(monosat.solver.Solve(a, b, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), result))
        self.assertFalse(monosat.solver.Solve(a, b, result))

    def test_nors(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        result = monosat.logic.Nor(a, b, c)
        self.assertFalse(monosat.solver.Solve(a, b, c, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), c, result))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(c), result))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c, result))
        self.assertFalse(monosat.solver.Solve(a, b, c, result))

    def test_xor(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        result = monosat.logic.Xor(a, b)
        self.assertFalse(monosat.solver.Solve(a, b, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(a, b, result))

    def test_xors(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        result = monosat.logic.Xor(a, b, c)
        self.assertTrue(monosat.solver.Solve(a, b, c, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), c, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(a, b, c, result))

    def test_xnor(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        result = monosat.logic.Xnor(a, b)
        self.assertTrue(monosat.solver.Solve(a, b, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(a, b, result))

    def test_xnors(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        result = monosat.logic.Xnor(a, b, c)
        self.assertFalse(monosat.solver.Solve(a, b, c, result))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), c, result))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c, monosat.logic.Not(result)))
        self.assertFalse(monosat.solver.Solve(a, b, c, result))

    def test_assertTrue(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        monosat.logic.AssertTrue(a)
        self.assertTrue(monosat.solver.Solve(a))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a)))

    def test_assertFalse(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        monosat.logic.AssertFalse(a)
        self.assertFalse(monosat.solver.Solve(a))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))

    def test_assertAnd(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        monosat.logic.AssertAnd(a, b)
        self.assertTrue(monosat.solver.Solve(a, b))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(b)))

    def test_assertAnds(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        monosat.logic.AssertAnd(a, b, c)
        self.assertTrue(monosat.solver.Solve(a, b, c))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(c)))

    def test_assertOr(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        monosat.logic.AssertOr(a, b)
        self.assertTrue(monosat.solver.Solve(a, b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(b)))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b)))

    def test_assertOrs(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        monosat.logic.AssertOr(a, b, c)
        self.assertTrue(monosat.solver.Solve(a, b, c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(c)))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(c)))

    def test_assertNand(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        monosat.logic.AssertNand(a, b)
        self.assertFalse(monosat.solver.Solve(a, b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(b)))

    def test_assertNands(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        monosat.logic.AssertNand(a, b, c)
        self.assertFalse(monosat.solver.Solve(a, b, c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(c)))

    def test_assertNor(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        monosat.logic.AssertNor(a, b)
        self.assertFalse(monosat.solver.Solve(a, b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertFalse(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(b)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b)))

    def test_assertNors(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        monosat.logic.AssertNor(a, b, c)
        self.assertFalse(monosat.solver.Solve(a, b, c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertFalse(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(c)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), monosat.logic.Not(c)))

    def test_assertXor(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        monosat.logic.AssertXor(a, b)
        self.assertFalse(monosat.solver.Solve(a, b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), b))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b)))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(b)))

    def test_assertXors(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        monosat.logic.AssertXor(a, b)
        self.assertFalse(monosat.solver.Solve(a, b, c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), b, c))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c))
        self.assertTrue(monosat.solver.Solve(a, monosat.logic.Not(b), c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(c)))

    def test_assertXnor(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        monosat.logic.AssertXnor(a, b)
        self.assertTrue(monosat.solver.Solve(a, b))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b)))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(b)))

    def test_assertXnors(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        monosat.logic.AssertXnor(a, b)
        self.assertTrue(monosat.solver.Solve(a, b, c))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), b, c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b), c))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b), c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(c))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(c)))

    def test_assertEqual(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        monosat.logic.AssertEq(a, b)
        self.assertTrue(monosat.solver.Solve(a, b))
        self.assertFalse(monosat.solver.Solve(monosat.logic.Not(a), b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a), monosat.logic.Not(b)))
        self.assertFalse(monosat.solver.Solve(a, monosat.logic.Not(b)))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(a)))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(monosat.logic.Not(b)))

    def test_assertAtMostOne(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        monosat.pbtheory.AssertAtMostOne([a, b, c])

        self.assertFalse(monosat.solver.Solve(a, b))
        self.assertFalse(monosat.solver.Solve(b, c))
        self.assertFalse(monosat.solver.Solve(a, c))
        self.assertTrue(monosat.solver.Solve(a, b.Not()))
        self.assertTrue(monosat.solver.Solve(b))
        self.assertTrue(monosat.solver.Solve(c), True)

    def test_unsatCore(self):
        monosat.Monosat().newSolver()
        a = monosat.logic.Var()
        b = monosat.logic.Var()
        c = monosat.logic.Var()
        d = monosat.logic.Var()
        x = monosat.logic.Equal(a, b)
        y = monosat.logic.Equal(b, c.Not())
        z = monosat.logic.Equal(c, a)
        w = monosat.logic.Equal(b, d.Not())
        q = monosat.logic.Equal(d, a)
        self.assertTrue(monosat.solver.Solve())
        self.assertFalse(monosat.solver.Solve(x, y, z, w, q))
        conflict = monosat.solver.getConflictClause()
        self.assertFalse(len(conflict) == 0)
        test = []
        for l in conflict:
            assert (not l.sign())
            test.append(l.Not())

        self.assertFalse(monosat.solver.Solve(test))

        self.assertFalse(monosat.solver.Solve(x, y, z, w, q))
        conflict2 = monosat.solver.getConflictClause(True)
        self.assertFalse(len(conflict2) == 0)
        self.assertTrue(len(conflict2) <= len(conflict))
        test2 = []
        for l in conflict2:
            assert (not l.sign())
            test2.append(l.Not())

        self.assertFalse(monosat.solver.Solve(test2))
        self.assertTrue(monosat.solver.Solve())

        conflict3 = monosat.solver.minimizeUnsatCore([x, y, z, w, q])
        self.assertFalse(len(conflict3) == 0)
        self.assertTrue(len(conflict3) <= 3)
        self.assertTrue(monosat.solver.Solve())


if __name__ == '__main__':
    unittest.main()
