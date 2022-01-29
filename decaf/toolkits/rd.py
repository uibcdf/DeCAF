# -*- coding: utf-8 -*-
"""RDKit toolkit for DeCAF"""

from decaf import PHARS, Pharmacophore
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from collections import deque
import math


PATTERNS = {phar: Chem.MolFromSmarts(smarts) for (phar, smarts) in PHARS.items()}


def __count_bonds(ligand, a1, a2, exclude):

    """Count number of bonds between two pharmacophore points, if the shortest
    path does not contain any other pharmacophore point.

    Args:
       a1, a2 (int): source and target atoms
       exclude (list): atoms (ids) that cannot be in the shortest path

    Returns:
       int: number of bonds in path or -1 if there is no path between a1 and a2
    """

    visited = []
    bonds_nr = -1
    queue = deque([(a1, 0)])
    while queue:
        idx, depth = queue.popleft()
        visited.append(idx)
        if idx == a2:
            bonds_nr = depth
            break
        else:
            for atom in ligand.GetAtomWithIdx(idx).GetNeighbors():
                if atom.GetIdx() not in visited and atom.GetIdx() not in exclude:
                    queue.append((atom.GetIdx(), depth+1))
    return bonds_nr


def phar_from_mol(ligand):
    """Create Pharmacophore from given RDKit.Chem.Mol object."""
    if not isinstance(ligand, Chem.Mol):
        raise TypeError("Invalid ligand! Expected RDKit.Chem.Mol object, got "
                        "%s instead" % type(ligand).__name__)
    matches = {}
    for (phar, pattern) in PATTERNS.items():
        atoms = list(zip(*ligand.GetSubstructMatches(pattern)))
        if len(atoms) > 0:
            matches[phar] = list(atoms[0])
        else:
            matches[phar] = []
    points = {}   # graph ids of matched atoms

    nodes = []
    idx = 0
    for (phar, atoms) in matches.items():
        for atom in atoms:
            if atom in points:
                nodes[points[atom]]["type"][phar] = 1.0
            else:
                nodes.append({"label": idx, "type": {phar: 1.0},
                             "freq": 1.0})
                points[atom] = idx
                idx += 1

    edges = np.zeros((idx, idx))

    keys = sorted(points.keys())
    for i in range(len(keys)):
        for j in range(i):
            dist = float(__count_bonds(ligand, keys[i], keys[j],
                         [keys[k] for k in range(len(keys)) if
                          k not in [i, j]]))
            if dist > -1:
                edges[points[keys[i]], points[keys[j]]] = dist
                edges[points[keys[j]], points[keys[i]]] = dist

    if not ligand.HasProp("_Name"):
        return Pharmacophore(nodes, edges, molecules=1.0)

    else:
        return Pharmacophore(nodes, edges, molecules=1.0,
                             title=ligand.GetProp("_Name"))


def layout(p):
    """Calculate points positions for depiction of Pharmacophore p using RDKit."""
    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore object, got %s instead" %
                        type(p).__name__)
    positions = np.zeros((p.numnodes, 2))
    m = Chem.MolFromSmarts("*")
    med = Chem.rdchem.EditableMol(m)
    for i in range(p.numnodes - 1):
        med.AddAtom(Chem.rdchem.Atom("*"))
    for i in range(p.numnodes):
        for j in range(i):
            if p.edges[i, j] > 0:
                tmp = int(math.ceil(p.edges[i, j])) - 1
                prev = i

                # add invisible atoms to get right distance
                for k in range(tmp):
                    idx = med.AddAtom(Chem.rdchem.Atom("*"))
                    med.AddBond(prev, idx)
                    prev = idx

                med.AddBond(prev, j)
    m = med.GetMol()
    for atom in m.GetAtoms():
        if atom.GetIdx() >= p.numnodes:
            atom.SetHybridization(Chem.rdchem.HybridizationType.SP)

    AllChem.Compute2DCoords(m, bondLength=1.0)
    c = m.GetConformer()

    for i in range(p.numnodes):
        positions[i][0] = c.GetAtomPosition(i).x
        positions[i][1] = c.GetAtomPosition(i).y
    return positions
