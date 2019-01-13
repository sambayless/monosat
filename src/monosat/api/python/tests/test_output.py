import unittest
import warnings

import monosat
import shutil, tempfile
import os


class TestOutput(unittest.TestCase):
    def test_output(self):
        # To write to an output file, create a new solver and pass it the name of a file to write to
        # All constraints created in that new solver will be copied to that file
        test_file = tempfile.NamedTemporaryFile()
        filename = test_file.name
        test_file.close()
        monosat.Monosat().newSolver(output_file=filename)
        a = monosat.Var()
        b = monosat.Var()
        c = monosat.Or(a, monosat.Not(b))

        monosat.Assert(c)
        monosat.Assert(monosat.Not(c))

        # trivial contradition, so result should be UNSAT
        result = monosat.Solve()
        self.assertFalse(result)

        monosat.Monosat().newSolver()  # create a new solver, this time without any output file set.
        # solver now has no constraints, so it should be SAT
        self.assertTrue(monosat.Solve())

        # read in the previously saved constraints (these can also be read in on the command line, eg ./monosat tutorial.gnf
        monosat.Monosat().loadConstraints(filename)
        result2 = monosat.Solve()
        self.assertEqual(result, result2)
        os.remove(filename)

    def test_setOutputFile(self):
        # this tests the deprecated setOutputFile interface
        test_file = tempfile.NamedTemporaryFile()
        filename = test_file.name
        test_file.close()
        monosat.Monosat().newSolver(arguments="-decide-theories")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            monosat.Monosat().setOutputFile(output_file=filename)
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            self.assertTrue("deprecated" in str(w[-1].message))

        a = monosat.Var()
        b = monosat.Var()
        c = monosat.Or(a, monosat.Not(b))

        monosat.Assert(c)
        monosat.Assert(monosat.Not(c))

        # trivial contradition, so result should be UNSAT
        result = monosat.Solve()
        self.assertFalse(result)

        monosat.Monosat().newSolver()  # create a new solver, this time without any output file set.
        # solver now has no constraints, so it should be SAT
        self.assertTrue(monosat.Solve())

        # read in the previously saved constraints (these can also be read in on the command line, eg ./monosat tutorial.gnf
        monosat.Monosat().loadConstraints(filename)
        result2 = monosat.Solve()
        self.assertEqual(result, result2)
        os.remove(filename)


if __name__ == '__main__':
    unittest.main()
