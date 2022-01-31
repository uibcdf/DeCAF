# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 10:11:53 2015

@author: martasd
"""
import sys
from pybel import readfile
from multiprocessing import Pool

target = sys.argv[1]
fp = sys.argv[2]

if len(sys.argv) > 3:
    cpus = int(sys.argv[3])
else:
    cpus = 12


IN_PATH = "../data/smi_files/"
OUT_PATH = "../data/results/" + fp + "/"

ref_file = IN_PATH + target + ".ism"
drugs_file = IN_PATH + "drugs.ism"


out_file = OUT_PATH + target + ".tab"

dataset = {}
dataset["ref"] = [(mol.title, mol.calcfp(fptype=fp))
                  for mol in readfile("smi", ref_file)]
dataset["drugs"] = [(mol.title, mol.calcfp(fptype=fp))
                    for mol in readfile("smi", drugs_file)]


def compute_scores(idx):
    name = dataset["drugs"][idx][0]
    res = name
    for p in dataset["ref"]:
        s = p[1] | dataset["drugs"][idx][1]
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
