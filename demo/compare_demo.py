# -*- coding: utf-8 -*-
"""
Compare active and inactive molecules to reference set.

Created on Mon Mar 16 10:11:56 2015
@author: Marta Stepniewska
"""

from pybel import readfile
from decaf.toolkits.ob import phar_from_mol
from decaf.utils import similarity
from multiprocessing import Pool
from time import sleep


ref = [mol for mol in readfile("smi", "ref.ism")]
act = [mol for mol in readfile("smi", "actives.ism")]
inact = [mol for mol in readfile("smi", "inactives.ism")]

prefix = "compare-demo"

dataset = {}
dataset["ref"] = [phar_from_mol(mol) for mol in ref]
dataset["act"] = [phar_from_mol(mol) for mol in act]
dataset["inact"] = [phar_from_mol(mol) for mol in inact]


def compute_scores(tup):
    name = dataset[tup[0]][tup[1]].title
    res = tup[0]+","+name
    for p in dataset["ref"]:
        res += ","+str(similarity(p, dataset[tup[0]][tup[1]])[0])
    return res


sets = []

print "Comparing %s active and %s inactive molecules to %s reference molecules." \
      % (len(act), len(inact), len(ref))


for moltype in ["act", "inact"]:
    for i in xrange(len(dataset[moltype])):
        sets.append((moltype, i))

pool = Pool()
results = pool.map_async(compute_scores, sets, chunksize=1)

num = float(len(act)+len(inact))
while not results.ready():
    print("{}% of comparisons left".format(round(100.0 * results._number_left / num, 2)))
    sleep(10)

results = results.get()

name = prefix+".csv"
print "Saving results to %s file" % name

f = open(name, "w")
f.write("label,name")
for i in xrange(len(ref)):
    f.write(",ref"+str(i))
f.write("\n")
for row in results:
    f.write(row+"\n")
f.close()
print "finished"
