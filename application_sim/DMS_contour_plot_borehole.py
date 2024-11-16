'''
TMI anomaly contour plotting script
===========================================

This script is used to plot contours along y = 0m.
All output files desired should be placed in the
same folder. This script should be used to plot 
the borehole simulation results.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append(r'C:\Users\Tomas\simpeg')
# input plot title here
plot_title = "Borehole model profile at y = 0 m" 

# files to work with
dir_path = os.getcwd() + '\\outputs\\sw_borehole\\'
file_list = os.listdir(dir_path)

# plot all results
fig = plt.figure(figsize=(10, 6))
plt.ylabel('TMI anomaly value (nT)',size=20)
plt.xlabel('x (m)',size=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title(plot_title,size=20)
plt.grid()

topo_filename = dir_path + "magnetics_topo" + ".txt"
data_filename = dir_path + "magnetics_data" + ".obs"

topo_xyz = np.loadtxt(str(topo_filename))
dobs = np.loadtxt(str(data_filename))

receiver_locations = dobs[:, 0:3]
dobs = dobs[:, -1]

x = receiver_locations[:,0]
y = receiver_locations[:,1]
ydata1 = []
ydata2 = []
for ii in range(0,len(y)):
    if y[ii] == 0.25:
        ydata1.append((x[ii],dobs[ii]))
    elif y[ii] == -0.25:
        ydata2.append((x[ii],dobs[ii]))

ydata = (np.array(ydata1) + np.array(ydata2))/2
xvals = [x for x, y in ydata]
cdata = [y for x, y in ydata]

plt.plot(xvals[1:-1],cdata[1:-1],marker='.',linestyle='-')
plt.show()