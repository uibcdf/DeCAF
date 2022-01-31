# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:15:00 2015

@author: martasd
"""

from random import sample, seed
from math import log, ceil
from pybel import Smarts, readstring

pept_bond = Smarts("[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]")


seed(2345)
min_size = 50
max_size = 10000
intervals = 100
rep = 100

source_file = "../data/chembl_20.ism"

f = open(source_file)
lines = f.readlines()
f.close()

n = len(lines)

log_min = log(min_size, 10)
log_max = log(max_size, 10)
interval_size = (log_max - log_min) / (intervals-1)


def bad(x):
    mol = readstring("smi", x)
    try:
        return (len(pept_bond.findall(mol)) > 20)
    except:
        return False


for i in xrange(intervals):
    sample_size = int(ceil(10**(log_min + i*interval_size)))
    sample_set = sample(lines, sample_size)
    f = open(("../data/statistical_model/samples/sample"+str(i)
              +"_"+str(sample_size)+".ism"), "w")
    for l in sample_set:
        while bad(l):
            l = sample(lines, 1)[0]
        f.write(l)
    f.close()


test_set = sample(lines, rep)
f = open(("../data/statistical_model/samples/test.ism"), "w")
for l in test_set:
    f.write(l)
f.close()
