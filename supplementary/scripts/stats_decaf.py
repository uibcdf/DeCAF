# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 10:11:53 2015

@author: martasd
"""
import sys
from decaf.utils import similarity
from decaf.toolkits import ob
from pybel import readfile
from multiprocessing import Pool


sample = sys.argv[1]

if len(sys.argv) > 2:
    cpus = int(sys.argv[2])
else:
    cpus = 12

IN_PATH = "../data/statistical_model/samples/"
OUT_PATH = "../data/statistical_model/results/DeCAF/"

ref_file = IN_PATH + sample + ".ism"
test_file = IN_PATH + "test.ism"


out_file = OUT_PATH + sample + ".tab"

dataset = {}
dataset["ref"] = [ob.phar_from_mol(mol) for mol in readfile("smi", ref_file)]
dataset["test"] = [ob.phar_from_mol(mol) for mol in readfile("smi", test_file)]


def compute_scores(idx):
    name = dataset["test"][idx].title
    res = name
    for p in dataset["ref"]:
        s = similarity(p, dataset["test"][idx], dist_tol=0,
                       coarse_grained=False)[0]
        res += "\t"+str(s)
    return res

pool = Pool(cpus)
results = pool.map_async(compute_scores, range(len(dataset["test"])))

results = results.get()
pool.close()
pool.join()

f = open(out_file, "w")
f.write("name")
for i in dataset["ref"]:
    f.write("\t"+i.title)
f.write("\n")
for row in results:
    f.write(row+"\n")
f.close()
