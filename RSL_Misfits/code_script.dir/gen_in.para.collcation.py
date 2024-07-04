#!/usr/bin/evn python
# This code is to generate the input file in.para.collocation for different viscosity structure girds

from os import path
import numpy as np
import pdb

# Generate viscosity grids values
lims1_upper_mantle = 19
lims2_upper_mantle = 21.5
inv_upper_mantle = 0.1 #0.5
visc_upper_mantle = np.arange(lims1_upper_mantle, lims2_upper_mantle+0.1, inv_upper_mantle)
print("upper_mantle_grids: {len(visc_upper_mantle)}")

# Generate viscosity grids values
lims1_lower_mantle = 20.5
lims2_lower_mantle = 23.5
inv_lower_mantle = 0.1 #0.5
visc_lower_mantle = np.arange(lims1_lower_mantle, lims2_lower_mantle+0.1, inv_lower_mantle)
print("lower_mantle_grids: {len(visc_lower_mantle)}")

# Generage nlog1 for l=1, l=2 and l=3 to 10
n_lower_mantle = len(visc_lower_mantle)

#pdb.set_trace()
# l1 : constant for each constant upper mantle visc, increase by inv_l1 for increase lower mantle visc
lim1_l1 = 6.0
lim2_l1 = 8.0
inv_l1 = (lim2_l1 - lim1_l1)/len(visc_upper_mantle) #0.5
val_l1 = np.arange(lim1_l1, lim2_l1, inv_l1)
val_l1_full = np.array([val_l1] * len(visc_lower_mantle)).flatten(order="F")

# l2 : increase by inv_l2 for each constant upper mantle visc, constant for each constant lower mantle
lim1_l2 = 4.5
lim2_l2 = 6.5
inv_l2 = (lim2_l2 - lim1_l2)/len(visc_lower_mantle)
val_l2 = np.arange(lim1_l2, lim2_l2, inv_l2)
val_l2_full = np.array([val_l2] * len(visc_upper_mantle)).flatten(order="C")

# l3 to 10
val_l3_full = np.array([7.0] * (len(visc_upper_mantle) * len(visc_lower_mantle))).flatten(order="C")  # set 7.0 for all the grids

ddir = "../fn_para_collocation.dir"
count = 0
for i_upper_mantle in range(len(visc_upper_mantle)):
    for i_lower_mantle in range(len(visc_lower_mantle)):
        fn = "in.para.collocation" + str(count+1).zfill(3)
        f = open(path.join(ddir, fn), 'w')
        f.write('%e\n'%10**visc_upper_mantle[i_upper_mantle])
        f.write('%e\n'%10**visc_lower_mantle[i_lower_mantle])
        f.write('%f\n'%val_l1_full[count])
        f.write('%f\n'%val_l2_full[count])
        f.write('%f\n'%val_l3_full[count])

        count += 1

f.close()
