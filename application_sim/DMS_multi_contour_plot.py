'''
TMI anomaly contour plotting script
===========================================

This script is used to plot contours along y = 0m.
All output files desired should be placed in the
same folder. The ending of a file should be a number
followed by 'm'. For example, magnetics_data_30m.obs.
This script also works for different concentration outputs
of the borehole metal model.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import re

sys.path.append(r'C:\Users\Tomas\simpeg')
# input titles here
legend_title = r'Fe$^{2+}$ concentration (mg/L)' 
plt_title = "Borehole model profile comparison at y = 0 m"
# input directory path here
dir_path = os.getcwd() + '\\outputs\\Fe_bh\\'

def break_string_at_underscore_and_dot(s):
    # Use regex to split the string at both _ and .
    return re.split(r'[_\.\(m)]', s)

def convert_to_float(s):
    num = float(s)
    if num.is_integer():
        return int(num)
    else:
        return num


# files to work with
file_list = os.listdir(dir_path)
# retrieve all file suffixes
suffix = []
for file in file_list:
    s = break_string_at_underscore_and_dot(file)
    suffix.append(s[3]) 
# retrieve unique values in suffix list and order them in ascending order
suffix = list(set(suffix))
float_list = [convert_to_float(s.replace(',','.')) for s in suffix]
float_list = sorted(float_list)
suffix = [str(f) for f in float_list]
# plot all results
fig = plt.figure(figsize=(10, 6))
plt.ylabel('TMI anomaly value (nT)',size=20)
plt.xlabel('x (m)',size=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title(plt_title,size=20)
plt.grid()
for e in suffix:
    e = e.replace('.',',')
    topo_filename = dir_path + "magnetics_topo_" + e + "m.txt"
    data_filename = dir_path + "magnetics_data_" + e + "m.obs"

    topo_xyz = np.loadtxt(str(topo_filename))
    dobs = np.loadtxt(str(data_filename))

    receiver_locations = dobs[:, 0:3]
    dobs = dobs[:, -1]

    x = receiver_locations[:,0]
    y = receiver_locations[:,1]
    ydata = []
    for ii in range(0,len(y)):
        if y[ii] == 0:
            ydata.append((x[ii],dobs[ii]))
    xvals = [x for x, y in ydata]
    cdata = [y for x, y in ydata]

    e = e.replace(',','.')
    plt.plot(xvals[1:-1],cdata[1:-1],marker='.',linestyle='-',label=e)
    

plt.legend(title=legend_title)
plt.show()