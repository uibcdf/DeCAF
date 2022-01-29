# -*- coding: utf-8 -*-
"""
Unit tests for DeCAF module

Created on Thu Apr  2 15:57:49 2015
@author: Marta Stepniewska
"""
import unittest
import numpy as np


class PharmacophoreTests(unittest.TestCase):

    def setUp(self):
        from decaf import Pharmacophore
        nodes = [{"label": 0, "freq": 2.0, "type": {"HH": 2.0, "AR": 2.0, "R": 2.0}},
                 {"label": 1, "freq": 2.0, "type": {"HH": 2.0, "AR": 2.0, "R": 2.0}},
                 {"label": 2, "freq": 2.0, "type": {"HH": 2.0, "AR": 2.0, "R": 2.0}},
                 {"label": 3, "freq": 2.0, "type": {"HH": 2.0, "AR": 2.0, "R": 2.0}},
                 {"label": 4, "freq": 2.0, "type": {"HH": 2.0, "AR": 2.0, "R": 2.0}},
                 {"label": 5, "freq": 2.0, "type": {"AR": 2.0, "R": 2.0}},
                 {"label": 6, "freq": 2.0, "type": {"HA": 2.0}},
                 {"label": 7, "freq": 2.0, "type": {"HH": 2.0}},
                 {"label": 8, "freq": 1.0, "type": {"HA": 1.0, "HD": 1.0}},
                 {"label": 9, "freq": 1.0, "type": {"HA": 1.0, "HD": 1.0}}]

        edges = np.array([[0., 1., 0., 0., 0., 1., 0., 0., 0., 0.],
                          [1., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
                          [0., 1., 0., 1., 0., 0., 0., 0., 0., 0.],
                          [0., 0., 1., 0., 1., 0., 0., 0., 0., 0.],
                          [0., 0., 0., 1., 0., 1., 0., 0., 0., 0.],
                          [1., 0., 0., 0., 1., 0., 1., 0., 0., 0.],
                          [0., 0., 0., 0., 0., 1., 0., 2., 0., 0.],
                          [0., 0., 0., 0., 0., 0., 2., 0., 1., 1.],
                          [0., 0., 0., 0., 0., 0., 0., 1., 0., 1.],
                          [0., 0., 0., 0., 0., 0., 0., 1., 1., 0.]])

        self.phar = Pharmacophore(nodes, edges, molecules=2,
                                  title="test")

    def tearDown(self):
        self.phar = None

    def testCreate(self):
        self.assertEqual(self.phar.numnodes, len(self.phar.nodes))
        for i in range(self.phar.numnodes):
            for j in range(i):
                self.assertEqual(self.phar.edges[i, j], self.phar.edges[j, i],
                                 msg=("Array is asymetric! %s!=%s for i=%s, j=%s" %
                                      (self.phar.edges[i, j],
                                       self.phar.edges[j, i], i, j)))

    def testIter(self):
        i = 0
        for node in self.phar:
            self.assertEqual(node, self.phar.nodes[i])
            i += 1

    def testAddNode(self):
        node = {"label": "CH3", "freq": 1.0, "type": {"HH": 2.0}}
        num = self.phar.numnodes + 1.0
        nodes = self.phar.nodes+[node]

        self.phar.add_node(node)

        self.assertEqual(num, self.phar.numnodes)
        self.assertEqual(num, len(self.phar.edges))
        self.assertEqual(num, len(self.phar.edges[0]))
        self.assertEqual(nodes, self.phar.nodes)

    def testDelNode(self):
        from random import randint
        idx = randint(0, self.phar.numnodes - 1)

        num = self.phar.numnodes - 1.0
        nodes = self.phar.nodes[:idx]+self.phar.nodes[idx+1:]

        self.phar.remove_node(idx)

        self.assertEqual(num, self.phar.numnodes)
        self.assertEqual(nodes, self.phar.nodes)

    def testAddEdge(self):
        l = 1.0
        num = np.sum(self.phar.edges > 0) / 2.0 + 1
        for idx1 in range(self.phar.numnodes):
            for idx2 in range(idx1):
                if self.phar.edges[idx1, idx2] == 0:
                    self.phar.add_edge(idx1, idx2, l)

                    self.assertEqual(num, np.sum(self.phar.edges > 0) / 2.0)
                    self.assertEqual(self.phar.edges[idx1, idx2],
                                     self.phar.edges[idx2, idx1])
                    self.assertEqual(self.phar.edges[idx1, idx2], l)
                    self.setUp()

    def testRemoveEdge(self):
        num = np.sum(self.phar.edges > 0) / 2.0 - 1.0
        for idx1 in range(self.phar.numnodes):
            for idx2 in range(idx1):
                if self.phar.edges[idx1, idx2] > 0:

                    self.phar.remove_edge(idx1, idx2)

                    self.assertEqual(self.phar.edges[idx1, idx2],
                                     self.phar.edges[idx2, idx1])
                    self.assertEqual(self.phar.edges[idx1, idx2], 0.0)
                    self.assertEqual(num, np.sum(self.phar.edges > 0) / 2.0)
                    self.setUp()

    def testSaveRead(self):
        from decaf import Pharmacophore
        from os import remove
        filename = "test.p"
        self.phar.save(filename)
        p_copy = Pharmacophore.read(filename)

        self.assertEqual(self.phar.numnodes, p_copy.numnodes)
        self.assertEqual(self.phar.nodes, p_copy.nodes)
        for i in range(p_copy.numnodes):
            for j in range(p_copy.numnodes):
                self.assertEqual(self.phar.edges[i, j], p_copy.edges[i, j])
        self.assertEqual(self.phar.title, p_copy.title)
        self.assertEqual(self.phar.molecules, p_copy.molecules)
        remove(filename)
        self.assertRaises(IOError, Pharmacophore.read, filename)
        self.assertRaises(IOError, Pharmacophore.save, p_copy, "nonexist/"+filename)

    def testValidation(self):
        from decaf import Pharmacophore

        self.assertRaises(TypeError, Pharmacophore, "a", self.phar.edges)
        self.assertRaises(TypeError, Pharmacophore, self.phar.nodes, "a")
        self.assertRaises(TypeError, Pharmacophore, self.phar.nodes,
                          self.phar.edges, molecules="a")
        self.assertRaises(ValueError, Pharmacophore, self.phar.nodes,
                          self.phar.edges, molecules=-1)
        self.assertRaises(TypeError, Pharmacophore, self.phar.nodes,
                          self.phar.edges, title=1)

        invalid = [([{"freq": 2.0, "type": {"HH": 2.0, "AR": 2.0}}] +
                    self.phar.nodes[1:], self.phar.edges),
                   ([{"label": 0, "type": {"HH": 2.0, "AR": 2.0}}] +
                    self.phar.nodes[1:], self.phar.edges),
                   ([{"label": 0, "freq": 2.0}]+self.phar.nodes[1:],
                    self.phar.edges),
                   ([{"label": 0, "freq": 2.0, "type": {"H": 2.0, "AR": 2.0}}] +
                    self.phar.nodes[1:], self.phar.edges),
                   (self.phar.nodes, self.phar.edges[:3][:, :3])]

        for args in invalid:
            self.assertRaises(ValueError, Pharmacophore, *args)

        self.assertRaises(TypeError, self.phar.add_node, "1")
        self.assertRaises(ValueError, self.phar.add_node, {})
        self.assertRaises(TypeError, self.phar.remove_node, "1")
        self.assertRaises(ValueError, self.phar.remove_node, -1)
        self.assertRaises(ValueError, self.phar.remove_node, self.phar.numnodes)
        self.assertRaises(TypeError, self.phar.add_edge, "0", 1, 2)
        self.assertRaises(TypeError, self.phar.add_edge, 0, "1", 2)
        self.assertRaises(TypeError, self.phar.add_edge, 0, 1, "2")
        self.assertRaises(ValueError, self.phar.add_edge, 0, 0, 2)
        self.assertRaises(ValueError, self.phar.add_edge, -1, 0, 2)
        self.assertRaises(ValueError, self.phar.add_edge, 0, self.phar.numnodes, 2)
        self.assertRaises(ValueError, self.phar.remove_edge, -1, 0)
        self.assertRaises(TypeError, self.phar.remove_edge, "0", 1)
        self.assertRaises(TypeError, self.phar.remove_edge, 0, "1")
        self.assertRaises(ValueError, self.phar.remove_edge, 0, self.phar.numnodes)
