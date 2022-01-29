# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 09:07:48 2016

@author: martas
"""
import unittest
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))


class ToolkitsTests(unittest.TestCase):

    def setUp(self):
        self.string = "Nc1ccc(C(O)O)cc1	mol1"
        self.numnodes = 9
        self.numedges = 10
        self.types = {"AR": 6, "HH": 5, "HA": 3, "HD": 3, "R": 6}

    def testCreateOb(self):
        from pybel import readstring
        import decaf.toolkits.ob as ob
        mol = readstring("smi", self.string)
        phar = ob.phar_from_mol(mol)
        self.assertEqual(phar.numnodes, self.numnodes)
        self.assertEqual(np.sum(phar.edges > 0) / 2.0, self.numedges)

        types = {t: 0 for t in self.types}
        for i in range(phar.numnodes):
            for t in list(phar.nodes[i]["type"].keys()):
                types[t] += 1
        self.assertEqual(types, self.types)

    def testValidationOb(self):
        import decaf.toolkits.ob as ob
        self.assertRaises(TypeError, ob.phar_from_mol, "c1ccccc1")
        self.assertRaises(TypeError, ob.phar_from_mol, 2)
        self.assertRaises(TypeError, ob.layout, "c1ccccc1")
        self.assertRaises(TypeError, ob.layout, 2)

    def testCreateRd(self):
        from rdkit.Chem import MolFromSmiles
        import decaf.toolkits.rd as rd
        molstring, name = self.string.split()
        mol = MolFromSmiles(molstring)
        mol.SetProp("_Name", name)
        phar = rd.phar_from_mol(mol)
        self.assertEqual(phar.numnodes, self.numnodes)
        self.assertEqual(np.sum(phar.edges > 0) / 2.0, self.numedges)

        types = {t: 0 for t in self.types}
        for i in range(phar.numnodes):
            for t in list(phar.nodes[i]["type"].keys()):
                types[t] += 1
        self.assertEqual(types, self.types)

    def testValidationRd(self):
        import decaf.toolkits.rd as rd
        self.assertRaises(TypeError, rd.phar_from_mol, "c1ccccc1")
        self.assertRaises(TypeError, rd.phar_from_mol, 2)
        self.assertRaises(TypeError, rd.layout, "c1ccccc1")
        self.assertRaises(TypeError, rd.layout, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
