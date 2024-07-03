#!/usr/bin/env python
from os import path
import numpy as np
import pdb

ddir='/Users/kkx_pp/Documents/Project_GIA/test_ice6g_kang/ice6g_data'
year = np.loadtxt(path.join(ddir,'years.reverse'))
steps1 = 480.
steps2 = 20.
steps3 = 10.

#pdb.set_trace()
inv1 = (122.-26.)/steps1
year_new1 = np.arange(122.,26.,-1*inv1)

inv2 = (26.-25.)/steps2
year_new2 = np.arange(26.,21.,-1*inv2)

inv3 = (21.-20.5)/steps3
year_new3 = np.arange(21.,0.,-1*inv3)

#write the year_new in txt and save
#pdb.set_trace()
year_new = np.concatenate([year_new1,year_new2,year_new3])

np.savetxt(path.join(ddir,'years.reverse.new'),year_new)

