"""
Mine drainage/Lithium model
===================================================

This model can be used to simulate lithium and mine 
drainage scenarios in-situ. The output is a TMI anomaly plot
inside of the borehole.
"""

#########################################################################
# Import Modules
# --------------
#

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys

sys.path.append(r'C:\Users\Tomas\simpeg')
from discretize import TensorMesh
from discretize.utils import mkvc, active_from_xyz
from SimPEG.utils import plot2Ddata
from SimPEG import maps
from SimPEG.potential_fields import magnetics

#############################################
# User Input 
# ----------
#
# Define code directions and certain variables here.
# For further edits, you may need to edit the code below this section.
#
write_output = True       # True to output the forward model results
output_location = 'Fe_bh'  # output folder within 'outputs'
gradient = False          # True if you want a gradient boundary
simple = True             # True for simple model, False for complex
boundary = 0              # define the boundary for the simple model, this is not used in the complex model                              
metal_concentration = 1000  # mg/L  
borehole_side_length = 0.1 # side length of the square cross section of the borehole in meters 
iron_ions = True          # True for iron ions, False for lithium ions
plot_sensors = True       # True to plot the sensor locations

# Constants
if iron_ions:
    amu = 55.845           # g/mol
    p2 = 24                
else:
    amu = 6.9410           # g/mol  
    p2 = 15

suffix = f'{metal_concentration}m'            # suffix for the output folder
#############################################
# Topography
# ----------
#
# Surface topography is defined as an (N, 3) numpy array. We create it here but
# topography could also be loaded from a file.
#

[x_topo, y_topo] = np.meshgrid(np.linspace(-50, 50, 100), np.linspace(-borehole_side_length/2, borehole_side_length/2, 20))
z_topo =  np.sqrt((x_topo)**2+(y_topo)**2)
x_topo, y_topo, z_topo = mkvc(x_topo), mkvc(y_topo), mkvc(z_topo)
xyz_topo = np.c_[x_topo, y_topo, z_topo]

#############################################
# Defining the Survey
# -------------------
#
# Here, we define survey that will be used for the simulation. Magnetic
# surveys are simple to create. The user only needs an (N, 3) array to define
# the xyz locations of the observation locations, the list of field components
# which are to be modeled and the properties of the Earth's field.
#

# Define the observation locations as an (N, 3) numpy array or load them.
xr = np.linspace(-50, 50, 50)
yr = np.linspace(-borehole_side_length/2, borehole_side_length/2, 3)
x, y = np.meshgrid(xr, yr)
x, y = mkvc(x.T), mkvc(y.T)
fun_interp = LinearNDInterpolator(np.c_[x_topo, y_topo], z_topo)
z = - borehole_side_length/2*np.ones(len(y))  # Sensor location in m below surface.
receiver_locations = np.c_[x, y, z]

# Define the component(s) of the field we want to simulate as a list of strings.
# Here we simulation total magnetic intensity data.
components = ["tmi"]

# Use the observation locations and components to define the receivers. To
# simulate data, the receivers must be defined as a list.
receiver_list = magnetics.receivers.Point(receiver_locations, components=components)

receiver_list = [receiver_list]

# Define the inducing field H0 = (intensity [nT], inclination [deg], declination [deg])
inclination = 0
declination = 90 
strength = 560000

source_field = magnetics.sources.UniformBackgroundField(
    receiver_list=receiver_list,
    amplitude=strength,
    inclination=inclination,
    declination=declination,
)

# Define the survey
survey = magnetics.survey.Survey(source_field)


#############################################
# Defining a Tensor Mesh
# ----------------------
#
# Here, we create the tensor mesh that will be used for the forward simulation.
#

dh = borehole_side_length/20
dx = 2
hx = [(dx, 50)]
hy = [(dh, 20)]
hz = [(dh, 20)]
mesh = TensorMesh([hx, hy, hz], "CCN")


#############################################
# Defining a Susceptibility Model
# -------------------------------
#
# Here, we create the model that will be used to predict magnetic data
# and the mapping from the model to the mesh. The model
# consists of a susceptible sphere in a less susceptible host.
#

# Define susceptibility values for each unit in SI
background_susceptibility = 1e-5
fresh_susceptibility = -1.65e-5  #-9e-6
concentration = metal_concentration/(amu*1000)   # Moles/L 
c = concentration*1000 # moles/m^3
cmol = 1.571e-6*p2
tt = 295 # temperature in K
iron_susceptibility =  fresh_susceptibility + c*cmol/tt  # iron measured: -1.575e-5
mixed_susceptibility1 = -1.625e-5
mixed_susceptibility2 = -1.6e-5

# Find cells that are active in the forward modeling (cells below surface)
ind_active = active_from_xyz(mesh, xyz_topo)


# Define the model
model = background_susceptibility * np.ones(ind_active.sum())

# Define mapping from model to active cells
nC = int(ind_active.sum())
model_map = maps.IdentityMap(nP=nC)  # model is a vlue for each active cell

# Define model blocks
ii = 0
if simple:
    xr = [boundary]
for e in xr: # for simple boundary with gradient, replace xr with [boundary] including the brackets

    ind_fresh = (
        (mesh.gridCC[ind_active, 0] >= min(x))
        & (mesh.gridCC[ind_active, 0] <= e)
    )
    if gradient:
        ind_mix1 = (
            (mesh.gridCC[ind_active, 0] > e)
            & (mesh.gridCC[ind_active, 0] <= e+10)
        )

        ind_mix2 = (
            (mesh.gridCC[ind_active, 0] > e+10)
            & (mesh.gridCC[ind_active, 0] <= e+30)
        )

        ind_iron= (
            (mesh.gridCC[ind_active, 0] > e+30)
            & (mesh.gridCC[ind_active, 0] <= max(x))
        )

        model[ind_mix1] = mixed_susceptibility1
        model[ind_mix2] = mixed_susceptibility2
    else:
        ind_iron = (
            (mesh.gridCC[ind_active, 0] > e)
            & (mesh.gridCC[ind_active, 0] <= max(x))
        )


    model[ind_iron] = iron_susceptibility
    model[ind_fresh] = fresh_susceptibility
    ii += -2

# Plot Model
fig = plt.figure(figsize=(9, 4))

plotting_map = maps.InjectActiveCells(mesh, ind_active, np.nan)
ax1 = fig.add_axes([0.1, 0.12, 0.73, 0.78])
mesh.plot_slice(
    plotting_map * model,
    normal="Y",
    ax=ax1,
    ind=int(mesh.shape_cells[1] / 2),
    grid=True,
    clim=(np.min(model), np.max(model))
)
ax1.set_title("Borehole model slice at y = 0 m")
ax1.set_xlabel("x (m)")
ax1.set_ylabel("z (m)")
if plot_sensors:
    p1 = ax1.plot(x,z,'o',label='sensors')
    ax1.legend(handles= p1)

ax2 = fig.add_axes([0.85, 0.12, 0.05, 0.78])
norm = mpl.colors.Normalize(vmin=np.min(model), vmax=np.max(model)) 
cbar = mpl.colorbar.ColorbarBase(ax2, norm=norm, orientation="vertical")
cbar.set_label("Magnetic Susceptibility (SI)", rotation=270, labelpad=15, size=12)

plt.show()


###################################################################
# Simulation: TMI Data for a Susceptibility Model
# -----------------------------------------------
#
# Here we demonstrate how to predict magnetic data for a magnetic
# susceptibility model using the integral formulation.
#

# Define the forward simulation. By setting the 'store_sensitivities' keyword
# argument to "forward_only", we simulate the data without storing the sensitivities
simulation = magnetics.simulation.Simulation3DIntegral(
    survey=survey,
    mesh=mesh,
    model_type="scalar",
    chiMap=model_map,
    ind_active=ind_active,
    store_sensitivities="forward_only"
)

# Compute predicted data for a susceptibility model
dpred = simulation.dpred(model)

# Plot
fig = plt.figure(figsize=(6, 5))
v_max = np.max(np.abs(dpred))

ax1 = fig.add_axes([0.1, 0.2, 0.8, 0.85])
plot2Ddata(
    receiver_list[0].locations,
    dpred,
    ax=ax1,
    ncontour=30,
    clim=(-v_max, v_max),
    contourOpts={"cmap": "bwr"},
)


ax1.set_title("Borehole TMI Anomaly")
ax1.set_xlabel("x (m)")
ax1.set_ylabel("y (m)")
ax1.tick_params(axis='y', labelsize=12)
ax1.grid(axis='y')
if plot_sensors:
    p1 = ax1.plot(x,y,'.',label='sensors')
    ax1.legend(handles= p1)

# y_axis components
ax1.yaxis.set_ticks(np.linspace(-borehole_side_length, borehole_side_length, 5))
ax1.set_aspect(aspect=80)

ax2 = fig.add_axes([0.15, 0.3, 0.7, 0.03])  # [left, bottom, width, height]
norm = mpl.colors.Normalize(vmin=-np.max(np.abs(dpred)), vmax=np.max(np.abs(dpred)))
cbar = mpl.colorbar.ColorbarBase(
    ax2, norm=norm, orientation='horizontal', cmap=mpl.cm.bwr
)
cbar.set_label("$nT$", rotation=0, labelpad=15, size=12)

plt.show()
#######################################################
# Optional: Export Data
# ---------------------
#
# Write the data and topography
#

if write_output:
    dir_path = os.path.dirname(__file__).split(os.path.sep)
    dir_path.extend(["outputs"])
    dir_path.extend([output_location])
    dir_path = os.path.sep.join(dir_path) + os.path.sep

    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    fname = dir_path + "magnetics_topo_" +suffix+".txt"
    np.savetxt(fname, np.c_[xyz_topo], fmt="%.4e")

    np.random.seed(211)
    maximum_anomaly = np.max(np.abs(dpred))
    noise = 0.02 * maximum_anomaly * np.random.randn(len(dpred))
    fname = dir_path + "magnetics_data_" +suffix+".obs"
    np.savetxt(fname, np.c_[receiver_locations, dpred + noise], fmt="%.4e")

