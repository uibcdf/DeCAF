# -*- coding: utf-8 -*-
"""OpenBabel toolkit for DeCAF"""

from decaf import PHARS, Pharmacophore
import pybel
import openbabel as ob
import numpy as np
from collections import deque
import math


PATTERNS = {phar: pybel.Smarts(smarts) for (phar, smarts) in PHARS.items()}


def __count_bonds(a1, a2, exclude):
    """Count number of bonds between two pharmacophore points, if the shortest
    path does not contain any other pharmacophore point.

    Args:
       a1, a2 (OBAtom): source and target atoms
       exclude (list): atoms (ids) that cannot be in the shortest path

    Returns:
       int: number of bonds in path or -1 if there is no path between a1 and a2
    """
    visited = []
    bonds_nr = -1
    queue = deque([(a1, 0)])
    while queue:
        atom, depth = queue.popleft()
        idx = atom.GetIdx()
        visited.append(idx)
        if atom == a2:
            bonds_nr = depth
            break
        else:
            for atom in ob.OBAtomAtomIter(atom):
                if atom.GetIdx() not in visited and atom.GetIdx() not in exclude:
                    queue.append((atom, depth+1))
    return bonds_nr


def phar_from_mol(ligand):
    """Create Pharmacophore from given pybel.Molecule object."""
    if not isinstance(ligand, pybel.Molecule):
        raise TypeError("Invalid ligand! Expected pybel.Molecule object, got "
                        "%s instead" % type(ligand).__name__)

    matches = {}
    for (phar, pattern) in PATTERNS.items():
        atoms = list(zip(*pattern.findall(ligand)))
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
                nodes.append({"label": atom, "type": {phar: 1.0},
                             "freq": 1.0})
                points[atom] = idx
                idx += 1

    edges = np.zeros((idx, idx))

    keys = sorted(points.keys())
    for i in range(len(keys)):
        for j in range(i):
            dist = float(__count_bonds(ligand.atoms[keys[i]-1].OBAtom,
                         ligand.atoms[keys[j]-1].OBAtom,
                         [keys[k] for k in range(len(keys)) if
                          k not in [i, j]]))
            if dist > -1:
                edges[points[keys[i]], points[keys[j]]] = dist
                edges[points[keys[j]], points[keys[i]]] = dist

    if ligand.title == "":
        return Pharmacophore(nodes, edges, molecules=1.0)

    else:
        return Pharmacophore(nodes, edges, molecules=1.0, title=ligand.title)


def layout(p):
    """Calculate points positions for depiction of Pharmacophore p using OpenBabel."""
    if not isinstance(p, Pharmacophore):
        raise TypeError("Expected Pharmacophore object, got %s instead" %
                        type(p).__name__)

    positions = np.zeros((p.numnodes, 2))
    m = pybel.Molecule(ob.OBMol())
    for i in range(p.numnodes):
        m.OBMol.NewAtom()
    idx = p.numnodes + 1
    for i in range(p.numnodes):
        for j in range(i):
            if p.edges[i, j] > 0:
                tmp = int(math.ceil(p.edges[i, j])) - 1
                prev = i + 1

                # add invisible atoms to get right distance
                for k in range(tmp):
                    atom = m.OBMol.NewAtom(idx)
                    atom.SetHyb(1)
                    m.OBMol.AddBond(prev, idx, 1)
                    prev = idx
                    idx += 1
                m.OBMol.AddBond(prev, j + 1, 1)
    m.draw(show=False, update=True)

    for i in range(p.numnodes):
        positions[i][0] = m.atoms[i].coords[0]
        positions[i][1] = m.atoms[i].coords[1]
    return positions
