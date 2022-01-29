# -*- coding: utf-8 -*-
"""
Create Pharmacophore models from set of molecules.

Created on Mon Mar 16 10:11:53 2015
@author: Marta Stepniewska
"""

from pybel import readfile
from decaf.toolkits.ob import phar_from_mol
from decaf.utils import similarity, combine_pharmacophores
import numpy as np
from multiprocessing import Pool, Manager
from scipy.cluster.hierarchy import average as upgma
from time import sleep
import sys


ligands = [mol for mol in readfile("smi", "test.ism")]
num_mols = len(ligands)
phars = [phar_from_mol(mol) for mol in ligands]
dist = np.zeros((num_mols*(num_mols - 1.0)/2.0))
cutoff = 0.7
prefix = "demo-model"


def compute_score(idx):
    return similarity(phars[idx[0]], phars[idx[1]])


manager = Manager()

print "Comparing %s molecules to each other" % num_mols

sets = []
for i in xrange(num_mols):
    for j in xrange(i+1, num_mols):
        sets.append((i, j))

pool = Pool()
res = pool.map_async(compute_score, sets, chunksize=1)

num = float(num_mols*(num_mols - 1.0)/2.0)
while not res.ready():
    print("{}% of comparisons left".format(round(100.0 * res._number_left / num, 2)))
    sleep(30)

sim = res.get()
pool.close()

for i in xrange(len(dist)):
    if sim[i][0] == 0.0:
            dist[i] = float("inf")
    else:
        dist[i] = 1.0 / sim[i][0]

print "Clustering molecules"
clustering = upgma(dist)

name = prefix+"-clustering.csv"
print "Clustering saved to %s file" % name

np.savetxt(name, clustering)
#clustering = np.genfromtxt(name)

singles = len(phars)
phars += [None]*len(clustering)
phars = manager.list(phars)
models = []


def merge((idx, row)):
    global phars, singles
    sys.stdout.flush()
    c = 1.0/row[2]
    i = int(row[0])
    j = int(row[1])
    if c >= cutoff:
        while phars[i] is None or phars[j] is None:
            sleep(0.1)
        phars[singles+idx] = combine_pharmacophores(phars[i], phars[j],
                                                    freq_cutoff=0.2)
    return 0


def get_models(idx):
    global models, singles
    if phars[idx] is None:
        i = int(clustering[idx-singles][0])
        j = int(clustering[idx-singles][1])
        get_models(i)
        get_models(j)
    else:
        models.append(phars[idx])

print "Merging highly similar models."
pool = Pool()
res = pool.map_async(merge, zip(range(len(clustering)), clustering), chunksize=1)

while not res.ready():
    print("{}% to merge left".format(round(100.0 * res._number_left / len(clustering), 2)))
    sleep(30)

res = res.get()
pool.close()


get_models(len(phars)-1)


print "Created %s model(s) from %s molecules." % (len(models), len(ligands))
for m in xrange(len(models)):
    print "\tmodel nr %s from %s molecule(s)" % (m, models[m].molecules)
    models[m].save((prefix+"-"+str(m)+".p"))

print "Models saved to %s-<nr>.p files" % prefix
