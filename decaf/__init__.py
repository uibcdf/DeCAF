# -*- coding: utf-8 -*-
"""
DeCAF - Discrimination, Comparison, Alignment tool for small molecules.

Created on Tue Feb 10 15:51:16 2015
by Marta Stepniewska
"""

import numpy as np
import warnings

warnings.simplefilter('always', UserWarning)

# SMARTS definition of pharmacophore points:
PHARS = {"HH": "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]",  # hydrophobic
         "AR": "[a]",  # aromatic
         "HA": "[!$([#1,#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]",   # acceptor
         "HD": "[!$([#6,H0,-,-2,-3]),$([!H0;#7,#8,#9])]",    # donor
         "R": "[r]"}    # ring


# colors used for pharmacophore depiction
COLORS = {"HH": "#FFFF00",  # hydrophobic
          "AR": "#FF9900",  # aromatic
          "HA": "#6666FF",    # acceptor
          "HD": "#CC0000",   # donor
          "R": "#AA44AA"}   # ring


class Pharmacophore(object):
    """Graph-like object representing pharmacophore.

    Args:
       nodes (list): list of dicts representing nodes
          each node has:
             * label (any hashable): node label, e.g. atom index
             * freq (float): node frequency
             * type (dict): dictionary describing its pharmacophoric
               properties, where keys are points types and values are numbers
               of molecules with this property

             e.g. common part of C2H5NH2 and N(C2H5)3 molecules would be coded as:
                [{"label": 0, "freq": 2.0, "type": {"HA": 2.0, "HD": 1.0}},
                {"label": 1, "freq": 2.0, "type": {"HH": 2.0}}]

       edges (numpy array): symmetrical array with zeros at diagonal, where
         edges[i, j] is edge length, expressed in number of bonds in
         molecule (0.0 means that i and j are not connected)

        molecules (float): number of molecules used to create model

        title (str): pharmacophore description

    Example:
       Pharmacophore created manually from C2H5NH2 and N(C2H5)3 molecules:

       >>> nodes = [{"label": "N", "freq": 2.0, "type": {"HA": 2.0, "HD": 1.0}},
       ... {"label": "CH3", "freq": 2.0, "type": {"HH": 2.0}}]
       >>> edges = np.array([[0,2],[2,0]])
       >>> p = Pharmacophore(nodes=nodes, edges=edges, molecules=2,
       ... title="C2H5NH2+N(C2H5)3")
    """

    def __init__(self, nodes, edges, molecules=1.0, title="Pharmacophore"):

        if not isinstance(nodes, list):
            raise TypeError("Invalid nodes list!")

        for i in range(len(nodes)):
            if not isinstance(nodes[i], dict):
                raise TypeError("Invalid node %s!" % i)
            elif not Pharmacophore.check_node(nodes[i]):
                raise ValueError("Invalid node %s!" % i)

        if not isinstance(edges, np.ndarray):
            raise TypeError("Invalid edges array!")

        elif not Pharmacophore.check_edges(edges):
            raise ValueError("Invalid edges array!")

        elif len(nodes) != len(edges):
            raise ValueError("Size of edge array does not match nodes number!")

        if isinstance(molecules, (float, int)):
            if molecules <= 0.0:
                raise ValueError("Molecules number must be greater than 0!")

        else:
            raise TypeError("Molecules number must be a float or an int, not %s!"
                            % type(molecules).__name__)

        if not isinstance(title, str):
            raise TypeError("Pharmacophore title must be a string, not %s!" %
                            type(title).__name__)

        self.nodes = [m.copy() for m in nodes]
        self.edges = np.array(edges, copy=True)
        self.numnodes = len(self.nodes)
        self.molecules = float(molecules)
        self.title = title

        if self.numnodes == 0:
            warnings.warn("Pharmacophore is empty!")

    def __iter__(self):
        return iter(self.nodes)

    @staticmethod
    def check_node(node):
        """Check if node is valid.

        Args:
           node (dict): node to check

        Returns:
           bool: True if node is valid, False otherwise.
        """
        has_type = False
        has_label = False
        has_freq = False
        if "freq" in node:
            if isinstance(node["freq"], float):
                has_freq = True
            elif isinstance(node["freq"], int):
                has_freq = True
                node["freq"] = float(node["freq"])

        if "label" in node:
            has_label = True

        if "type" in node:
            has_type = True
            for key in node["type"]:
                if key not in PHARS:
                    has_type = False
                    break
                if isinstance(node["type"][key], int):
                    node["type"][key] = float(node["type"][key])
                elif not isinstance(node["type"][key], float):
                    has_type = False
                    break
        if has_type and has_label and has_freq:
            return True
        else:
            return False

    @staticmethod
    def check_edges(edges):
        """Check if edges array is valid.

        Args:
           edges (numpy array): array representing edges in the graph

        Returns:
           bool: True if array is valid, False otherwise.
        """
        if len(edges) == 0:
            return True

        if (edges.T == edges).all() and (edges.diagonal() == 0).all() and \
           np.min(edges) >= 0:
            return True
        else:
            return False

    def copy(self):
        """Generate deep copy of a Pharmacophore."""
        from copy import deepcopy
        return deepcopy(self)

    def add_edge(self, i, j, dist):
        """Add edge of length dist between nodes i and j.

        Args:
           i, j (int): nodes indices
           dist (float): edges length
        """
        if not isinstance(i, int) or not isinstance(j, int):
            raise TypeError("Node indicies must be int!")
        if not isinstance(dist, (int, float)):
            raise TypeError("Distance must be float or int!")
        if dist < 0:
            raise ValueError("Distance must be greater than or equal 0")
        elif i == j:
            raise ValueError("Cant add edge from node to itself!")
        elif i < 0 or i >= self.numnodes:
            raise ValueError("Begining node index out of range!")
        elif j < 0 or j >= self.numnodes:
            raise ValueError("End node index out of range!")
        else:
            self.edges[i, j] = self.edges[j, i] = dist

    def add_node(self, node):
        """Add node to Pharmacophore.

        Args:
           node (dict): node representation
        """
        if not isinstance(node, dict):
            raise TypeError("Invalid node!")
        if Pharmacophore.check_node(node):
            self.nodes.append(node.copy())
            self.edges = np.append(np.append(self.edges,
                                             np.zeros((1, self.numnodes)),
                                             axis=0),
                                   np.zeros((self.numnodes+1, 1)), axis=1)
            self.numnodes += 1
        else:
            raise ValueError("Invalid node!")

    def remove_node(self, i):
        """Remove node from model.

        Args:
           i (int): node index
        """
        if not isinstance(i, int):
            raise TypeError("Node id must be int!")
        if i < 0 or i >= self.numnodes:
            raise ValueError("Node index out of range!")
        else:
            del self.nodes[i]
            self.edges = np.delete(np.delete(self.edges, i, 0), i, 1)
            self.numnodes -= 1

            if self.numnodes == 0:
                warnings.warn("Last node removed. Pharmacophore is empty!")

    def remove_edge(self, i, j):
        """Remove edge between nodes i and j.

        Args:
           i, j (int): nodes indices
        """
        if not isinstance(i, int) or not isinstance(j, int):
            raise TypeError("Node indicies must be int!")
        if i < 0 or i >= self.numnodes:
            raise ValueError("Begining node index out of range!")
        elif j < 0 or j >= self.numnodes:
            raise ValueError("End node index out of range!")
        else:
            self.edges[i, j] = self.edges[j, i] = 0.0

    def save(self, filename):
        """Save Pharmacophore to a file using pickle module.

        Args:
           filename (str): path to a file
        """
        from pickle import dump
        try:
            with open(filename, "wb") as f:
                dump(self, f)
        except IOError as e:
            raise IOError("Cant write file!", e)

    @staticmethod
    def read(filename):
        """Read pickled Pharmacophore from a file.

        Args:
           filename (str): path to a file

        Returns:
           Pharmacophore: object read from a file
        """
        from pickle import load
        try:
            with open(filename, "rb") as f:
                return load(f)
        except IOError as e:
            raise IOError("Cant read file!", e)
