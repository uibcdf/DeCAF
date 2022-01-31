# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 09:41:51 2016

@author: martasd
"""

import sys
import numpy as np
from pybel import readfile
from openbabel import OBConformerSearch
from usrcat.toolkits.ob import generate_moments
from usrcat.sim import similarity
from multiprocessing import Pool

target = sys.argv[1]
mode = sys.argv[2]

if len(sys.argv) > 3:
    cpus = int(sys.argv[3])
else:
    cpus = 12


IN_PATH = "../data/smi_files/"
OUT_PATH = "../data/results/USRCAT/" + mode + "/"

ref_file = IN_PATH + target + ".ism"
drugs_file = IN_PATH + "drugs.ism"


out_file = OUT_PATH + target + ".tab"


mols = {}
for key, fname in {"ref": ref_file, "drugs": drugs_file}.iteritems():
    mols[key] = [mol for mol in readfile("smi", fname)]


def gen_moments((key, idx, mode)):
    mol = mols[key][idx]
    mol.make3D()
    if mode == 'localopt':
        mol.localopt()
        return (mol.title, generate_moments(mol))

    elif mode == 'all':
        obconf = OBConformerSearch()
        obconf.Setup(mol.OBMol)
        obconf.Search()
        obconf.GetConformers(mol.OBMol)
        mol.OBMol.NumConformers()
        mol_moments = []
        for i in xrange(mol.OBMol.NumConformers()):
            mol.OBMol.SetConformer(i)
            mol_moments.append(generate_moments(mol))
        return (mol.title, np.array(mol_moments))


pool = Pool(cpus)

dataset = {}
for key in mols:
    n = len(mols[key])
    results = pool.map_async(gen_moments, zip([key]*n, range(n), [mode]*n))
    dataset[key] = results.get()


def compute_scores(idx):
    name = dataset["drugs"][idx][0]
    res = name
    for p in dataset["ref"]:
        s = similarity(p[1], dataset["drugs"][idx][1])[1]
        res += "\t"+str(s)
    return res

pool = Pool(cpus)
results = pool.map_async(compute_scores, range(len(dataset["drugs"])))

results = results.get()
pool.close()
pool.join()

f = open(out_file, "w")
f.write("name")
for i in dataset["ref"]:
    f.write("\t"+i[0])
f.write("\n")
for row in results:
    f.write(row+"\n")
f.close()
