# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 07:07:28 2016

@author: martasd
"""

import sys
import pandas as pd
from glob import glob

path_pattern = sys.argv[1]

for fname in glob(path_pattern):
    table = pd.read_table(fname)
    compname = fname+'.bz2'
    table.to_csv(compname, index=False, sep='\t', compression='bz2')
