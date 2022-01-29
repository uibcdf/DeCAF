# -*- coding: utf-8 -*-
"""
Filter database using Pharmacophore model.

Created on Mon Mar 16 10:11:53 2015
@author: Marta Stepniewska
"""

from pybel import readfile
from decaf import Pharmacophore
from decaf.toolkits.ob import phar_from_mol
from decaf.utils import similarity
from multiprocessing import Process, Manager, cpu_count
from time import sleep


NUM_PROCESSES = cpu_count()
cutoff = 0.8


database = readfile("smi", "all.ism")
model = Pharmacophore.read("model.p")

print "Read model with %s nodes created from %s molecules." % (model.numnodes,
                                                               model.molecules)


manager = Manager()
similar = manager.list()
proc = [None]*NUM_PROCESSES


def check_mol(mol):
    p = phar_from_mol(mol)
    s, c = similarity(model, p)
    #print s, c
    if s > cutoff:
        similar.append((mol.write(), s, c))


i = 0
compared = 0
while True:
    if (proc[i] is not None and not proc[i].is_alive()) or proc[i] is None:
        try:
            mol = database.next()
            compared += 1
            if compared % 10 == 0:
                print compared, "molecules comapred"
            proc[i] = Process(target=check_mol, args=(mol,))
            proc[i].start()
        except StopIteration:
            break
    i = (i + 1) % NUM_PROCESSES

print "All", compared, "molecules comapred."

for p in proc:
    while p.is_alive():
        sleep(0.1)

print "Found %s similar molecules:" % len(similar)
for s in similar:
    print "SMILES:", s[0].split()[0]
    print "description:", s[0].split("\t")[1].rstrip()
    print "score: %s, cost: %s" % (s[1:])
    print ""
