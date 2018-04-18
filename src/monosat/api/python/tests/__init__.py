import unittest
import tests

def load_tests(loader, testset, ignore):
    testset.addTests(tests.TestLogic)
    testset.addTests(tests.TestGraph)
    testset.addTests(tests.TestOutput)
    return testset