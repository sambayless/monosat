import unittest
import monosat


class TestGraph(unittest.TestCase):

    def test_nNodes(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(10):
            g.addNode()
            self.assertEqual(g.nNodes(), i + 1)

    def test_nEdges(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with one diagonal edges
        #
        # 0 *--* 1
        #  |/ |
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1)
        e_0_2 = g.addEdge(0, 2)
        e_1_3 = g.addEdge(1, 3)
        e_1_2 = g.addEdge(1, 2)
        e_2_3 = g.addEdge(2, 3)

        self.assertEqual(g.nEdges(), 5)

    def test_addNode(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(10):
            g.addNode()
            self.assertEqual(g.nNodes(), i + 1)

    def test_reaches(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with one diagonal edge
        #
        # 0 *--* 1
        #  |/ |
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1)
        e_0_2 = g.addEdge(0, 2)
        e_1_3 = g.addEdge(1, 3)
        e_1_2 = g.addEdge(1, 2)
        e_2_3 = g.addEdge(2, 3)

        r = g.reaches(0, 3)
        r2 = g.reaches(0, 3)
        self.assertTrue(monosat.Solve(r))
        self.assertFalse(monosat.Solve(r, e_0_1.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))

        self.assertFalse(monosat.Solve(r, e_0_2.Not(), e_1_3.Not(), e_2_3.Not()))

        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))
        # There should only be one solution to this: 0->1, 1->2, 2->3
        nodes = g.getPath(r, False)
        edges = g.getPath(r, True)
        self.assertEqual(len(edges), 3)
        self.assertEqual(len(nodes), 4)

        self.assertEqual(edges[0], e_0_1)
        self.assertEqual(edges[1], e_1_2)
        self.assertEqual(edges[2], e_2_3)
        self.assertTrue(e_0_1.value())
        self.assertTrue(e_1_2.value())
        self.assertTrue(e_2_3.value())
        self.assertFalse(e_0_2.value())
        self.assertFalse(e_1_3.value())

        self.assertTrue(monosat.Solve(r))
        self.assertTrue(r.value())
        self.assertTrue(r2.value())
        self.assertEqual(g.getPath(r, False), g.getPath(r2, False))
        self.assertEqual(g.getPath(r, True), g.getPath(r2, True))

    def test_reachesBack(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with one diagonal edge
        #
        # 0 *--* 1
        #  |/ |
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1)
        e_0_2 = g.addEdge(0, 2)
        e_1_3 = g.addEdge(1, 3)
        e_1_2 = g.addEdge(1, 2)
        e_2_3 = g.addEdge(2, 3)

        r = g.reachesBackward(3, 0)
        r2 = g.reachesBackward(0, 3)
        self.assertTrue(monosat.Solve(r))
        self.assertFalse(monosat.Solve(r2))
        self.assertTrue(monosat.Solve(r))
        self.assertFalse(monosat.Solve(r, e_0_1.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))

        self.assertFalse(monosat.Solve(r, e_0_2.Not(), e_1_3.Not(), e_2_3.Not()))

        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))
        # There should only be one solution to this: 0->1, 1->2, 2->3
        nodes = g.getPath(r, False)
        edges = g.getPath(r, True)
        self.assertEqual(len(edges), 3)
        self.assertEqual(len(nodes), 4)

        self.assertEqual(edges[2], e_0_1)
        self.assertEqual(edges[1], e_1_2)
        self.assertEqual(edges[0], e_2_3)
        self.assertTrue(e_0_1.value())
        self.assertTrue(e_1_2.value())
        self.assertTrue(e_2_3.value())
        self.assertFalse(e_0_2.value())
        self.assertFalse(e_1_3.value())

        self.assertTrue(monosat.Solve(r))

    def test_onPath(self):
        monosat.Monosat().newSolver(output_file="/tmp/test.gnf")
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with one diagonal edge
        #
        # 0 *--* 1
        #  |/ |
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1)
        e_0_2 = g.addEdge(0, 2)
        e_1_3 = g.addEdge(1, 3)
        e_1_2 = g.addEdge(1, 2)
        e_2_3 = g.addEdge(2, 3)

        r = g.onPath(1, 0, 3)
        r2 = g.onPath(1, 0, 3)
        self.assertFalse(monosat.Solve(r, ~e_0_1))
        self.assertFalse(monosat.Solve(r, e_0_1.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))

        self.assertFalse(monosat.Solve(r, e_1_3.Not(), e_2_3.Not()))
        self.assertFalse(monosat.Solve(r, e_0_2.Not(), e_1_3.Not(), e_2_3.Not()))

        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))
        # There should only be one solution to this: 0->1, 1->2, 2->3
        nodes = g.getPath(r, False)
        edges = g.getPath(r, True)
        self.assertEqual(len(edges), 3)
        self.assertEqual(len(nodes), 4)

        self.assertEqual(edges[0], e_0_1)
        self.assertEqual(edges[1], e_1_2)
        self.assertEqual(edges[2], e_2_3)
        self.assertTrue(e_0_1.value())
        self.assertTrue(e_1_2.value())
        self.assertTrue(e_2_3.value())
        self.assertFalse(e_0_2.value())
        self.assertFalse(e_1_3.value())

        self.assertTrue(monosat.Solve(r))
        self.assertTrue(r.value())
        self.assertTrue(r2.value())
        self.assertEqual(g.getPath(r, False), g.getPath(r2, False))
        self.assertEqual(g.getPath(r, True), g.getPath(r2, True))

    def test_maximumFlow_geq(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with one diagonal edge
        #
        # 0 *--* 1
        #  |/ |
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1, monosat.BitVector(4, 1))
        e_0_2 = g.addEdge(0, 2, monosat.BitVector(4, 1))
        e_1_3 = g.addEdge(1, 3, monosat.BitVector(4, 1))
        e_1_2 = g.addEdge(1, 2, monosat.BitVector(4, 2))
        e_2_3 = g.addEdge(2, 3, monosat.BitVector(4, 1))
        cmp = monosat.BitVector(4)
        f = g.maxFlowGreaterOrEqualTo(0, 3, cmp)
        self.assertTrue(monosat.Solve(f))
        self.assertFalse(monosat.Solve(f, cmp.eq(3)))
        self.assertTrue(monosat.Solve(f, cmp.eq(2)))
        self.assertTrue(monosat.Solve(f, cmp.eq(1)))
        self.assertTrue(monosat.Solve(f, e_0_1.Not(), e_2_3.Not()))
        self.assertFalse(monosat.Solve(f, e_0_1.Not(), e_2_3.Not(), cmp.eq(1)))
        self.assertTrue(monosat.Solve(f, e_0_2.Not(), e_1_3.Not()))
        self.assertTrue(monosat.Solve(f))

    def test_distance(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with one diagonal edge
        #
        # 0 *--* 1
        #  |/ |
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1, monosat.BitVector(4, 1))
        e_0_2 = g.addEdge(0, 2, monosat.BitVector(4, 1))
        e_1_3 = g.addEdge(1, 3, monosat.BitVector(4, 1))
        e_1_2 = g.addEdge(1, 2, monosat.BitVector(4, 2))
        e_2_3 = g.addEdge(2, 3, monosat.BitVector(4, 1))

        dist = monosat.BitVector(4)
        d = g.distance_leq(0, 3, dist)
        monosat.AssertTrue(d)
        self.assertTrue(monosat.Solve(dist.gt(0)))
        self.assertFalse(monosat.Solve(dist.eq(0)))
        self.assertFalse(monosat.Solve(dist.eq(1)))
        self.assertTrue(monosat.Solve(dist.eq(2)))
        self.assertTrue(monosat.Solve(dist.eq(3)))

        self.assertFalse(monosat.Solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(1)))
        self.assertTrue(monosat.Solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(2)))
        self.assertTrue(monosat.Solve(e_0_1, e_0_2, e_1_3, e_2_3, dist.eq(3)))

        self.assertFalse(monosat.Solve(e_0_1.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(e_0_2.Not(), e_1_3.Not()))
        self.assertFalse(
            monosat.Solve(e_0_2.Not(), e_1_3.Not(), e_2_3.Not(), dist.eq(1))
        )

    def test_acyclicDirected(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with two diagonal edges
        #
        # 0 *--* 1
        #  |/\|
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1)
        e_0_2 = g.addEdge(0, 2)
        e_1_3 = g.addEdge(1, 3)
        e_1_2 = g.addEdge(1, 2)
        e_2_3 = g.addEdge(2, 3)
        e_3_0 = g.addEdge(3, 0)

        r = g.acyclic()
        self.assertTrue(monosat.Solve(r))
        self.assertTrue(monosat.Solve(r, e_0_1.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))
        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not(), e_2_3.Not()))
        self.assertFalse(monosat.Solve(r, e_0_1, e_1_2, e_2_3, e_3_0))
        self.assertTrue(monosat.Solve(r))

        self.assertTrue(monosat.Solve(r.Not()))
        self.assertFalse(monosat.Solve(r.Not(), e_0_1.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r.Not(), e_0_2.Not(), e_1_3.Not()))
        self.assertFalse(monosat.Solve(r.Not(), e_0_2.Not(), e_1_3.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r.Not(), e_0_1, e_1_2, e_2_3, e_3_0))
        self.assertTrue(monosat.Solve(r.Not()))
        self.assertTrue(monosat.Solve(r))

    def test_acyclicUndirected(self):
        monosat.Monosat().newSolver()
        g = monosat.Graph()
        for i in range(4):
            g.addNode()

        # create a directed square graph with two diagonal edges
        #
        # 0 *--* 1
        #  |/\|
        # 2 *--* 3

        e_0_1 = g.addEdge(0, 1)
        e_0_2 = g.addEdge(0, 2)
        e_1_3 = g.addEdge(1, 3)
        e_1_2 = g.addEdge(1, 2)
        e_2_3 = g.addEdge(2, 3)
        e_3_0 = g.addEdge(3, 0)

        r = g.acyclic(False)
        self.assertTrue(monosat.Solve(r))
        self.assertTrue(monosat.Solve(r, e_0_1.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not()))
        self.assertTrue(monosat.Solve(r, e_0_2.Not(), e_1_3.Not(), e_2_3.Not()))
        self.assertFalse(monosat.Solve(r, e_0_1, e_1_2, e_2_3, e_3_0))
        self.assertTrue(monosat.Solve(r))

        self.assertTrue(monosat.Solve(r.Not()))
        # valid undirected cycle: e_3_0 -> e_0_2 -> e_1_2 -> e_1_3
        self.assertTrue(monosat.Solve(r.Not(), e_0_1.Not(), e_2_3.Not()))
        self.assertFalse(monosat.Solve(r.Not(), e_0_1.Not(), e_2_3.Not(), e_1_3.Not()))
        self.assertTrue(monosat.Solve(r.Not(), e_0_2.Not(), e_1_3.Not()))
        self.assertFalse(monosat.Solve(r.Not(), e_0_2.Not(), e_1_3.Not(), e_2_3.Not()))
        self.assertTrue(monosat.Solve(r.Not(), e_0_1, e_1_2, e_2_3, e_3_0))
        self.assertTrue(monosat.Solve(r.Not()))
        self.assertTrue(monosat.Solve(r))


if __name__ == "__main__":
    unittest.main()
