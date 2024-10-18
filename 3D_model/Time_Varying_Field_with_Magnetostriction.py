# This code is used to run a single simulation with
# a 3D magnetostriction model.

# Run Number
num = 0             # define the number label of the run here, this will create a new folder for the outputs

# Import packages
import oommfc as mc
import discretisedfield as df
import micromagneticmodel as mm
import numpy as np
import matplotlib.pyplot as plt
import random
import os

# Time setup
tfinal = 1          # Final time of simulation
tstep = 1e-6        # timestep in seconds
ntot = 1000         # total number of time steps
# Material properties
A_ex = 15e-12       # exchange constant
L1 = -24            # Strain constant 1 in ppm
L2 = -48            # Strain constant 2 in ppm
#Ku = -0.5e4         # anisotropy constant
#u_easy = (1, 0, 1)  # easy axis
#u_hard = (0, 1, 0)  # hard axis
damping = 0.9       # damping constant (alpha)
Ms = 4.8e5          # magnetization saturation (A/m)
# Alternating source properties
source_freq = 100   # source frequency (Hz)
H_app = (0, 0, 1.3e4)# applied field (A/m)
source_type = 'sin' # waveform of source
# Strain measurement axis
B = np.array([0,0,1])
# Geometry of object
r = 5e-9            # material radius, m
L = 1e-8            # material length, m
d = .5e-9            # domain dimensions, m

# Create a new folder with name model_{num}
newpath = r'C:\Users\Tomas\Desktop\Martin Group\ubermag_images\model_'
newpath = newpath + str(num)
if not os.path.exists(newpath):
    os.makedirs(newpath)

# Create a new file with all parameters written into it
newfile = newpath + r'\parameters.txt'
f = open(newfile,"w")
f.write(
f"""
# Run Number
num = {str(num)}                   # define the number label of the run here, this will create a new folder for the outputs
# Time setup
tfinal = {str(tfinal)}             # Final time of simulation in seconds
tstep = {str(tstep)}               # timestep in seconds
ntot = {str(ntot)}                 # total number of time steps
# Material properties
A_ex = {str(A_ex)}                 # exchange constant
L1 = {str(L1)}                     # Strain constant 1 in ppm
L2 = {str{L2}}                     # Strain constant 2 in ppm
damping = {str(damping)}           # damping constant (alpha)
Ms = {str(Ms)}                     # magnetization saturation (A/m)
# Alternating source properties
source_freq = {str(source_freq)}   # source frequency (Hz)
H_app = {str(H_app)}               # applied field (A/m)
source_type = {source_type}        # waveform of source
# Strain measurement axis
B = {str(B)}
# Geometry of object
r = {str(r)}                       # material radius (m)
L = {str(L)}                       # material length (m)
d = {str(d)}                       # domain dimensions (m)
"""
)

f.close()

# Define initial mesh
p1 = (-r, -r, -L/2)                        # Starting point
p2 = (r, r, L/2)                           # Ending point
cell = (d, d, d)                           # Cell size
region = df.Region(p1=p1, p2=p2)           # Define the region
mesh = df.Mesh(region=region, cell=cell)   # Create the mesh

# Define the system name
system = mm.System(name='time_dependent_field')

# Define system energy
system.energy = (mm.Zeeman(H=H_app, func='sin', f=source_freq, t0=0)
                + mm.Demag()
                + mm.Exchange(A=A_ex))
#                + mm.CubicAnisotropy(K=Ku, u1=u_easy, u2=u_hard))

# Define system dynamics
system.dynamics = mm.Precession(gamma0=mm.consts.gamma0) + mm.Damping(alpha=damping)

def Ms_fun(pos):
    """Function to set magnitude of magnetisation: zero outside cylindric shape,
    Ms inside cylinder.

    Cylinder radius is r.

    """
    x, y, z = pos
    if (x**2 + y**2)**0.5 < r:
        return Ms
    else:
        return 0

# Define the field element
system.m = df.Field(mesh, dim=3, value=(0,0,0), norm=Ms_fun)

# Assign a random value of 4 defined states to each cell to represent the different domains.
i_states = {'state_1':(1,1,1),'state_2':(-1,-1,-1),'state_3':(-1,1,-1),'state_4':(1,-1,1)}
n = len(system.m.array[0,0])
i=0
for array_mass in system.m.array:
    j=0
    for array_element in system.m.array[i]:
        array_element = np.array(random.choices(list(i_states.values()),k=n))
        system.m.array[i,j] = array_element
        j+=1
    i+=1

# Redefine the field element to define the cylindrical shape again
system.m = df.Field(mesh, dim=3, value = system.m.array, norm=Ms_fun)

# Plot the initial magnetization of one cross section that is one cell thick
system.m.plane('z').mpl()
plt.savefig(newpath+r'\initial_mag.jpg')

# Evolve the system
ev = mc.RungeKuttaEvolver(min_timestep=tstep) # Define the evolver as Runge Kutta method
td = mc.TimeDriver(evolver=ev)                # Setup time driver with evolver input
td.drive(system, t=tfinal, n=ntot, verbose=2) # Drive the system

# Plot the external magnetic field
system.table.mpl(y=['Bx_zeeman', 'By_zeeman', 'Bz_zeeman'])
plt.legend(['$H_x$','$H_y$','$H_z$'])
plt.title('External Field')
plt.ylabel('Field Magnitude (kA/m)')
plt.savefig(newpath+r'\mag_field.jpg')

# Plot the normalized magnetization
system.table.mpl(y=['mx', 'my', 'mz'])
plt.legend(['$m_x$','$m_y$','$m_z$'])
plt.ylabel('m')
plt.title('Normalized Magnetization')
plt.savefig(newpath+r'\magnetization.jpg')

# Assign magnetization data to x,y,z variables.
mx = system.table.data['mx']
my = system.table.data['my']
mz = system.table.data['mz']

# Normalize the strain measurement axis and separate into x,y,z.
B = B/np.linalg.norm(B)
Bx = B[0]
By = B[1]
Bz = B[2]

# Strain calculation for cubic anisotropy
ll = 3/2*L1*(mx**2*Bx**2+my**2*By**2+mz**2*Bz**2-1/3)+3*L2*(mx*my*Bx*By+my*mz*By*Bz+mz*mx*Bz*Bx)

# Define a new column in the data for the strain
system.table.data['ll'] = ll

# Plot the strain response
system.table.mpl(y=['ll'])
plt.legend().remove()
plt.title('Model Strain Response')
plt.ylabel(r'$\lambda$ (ppm)')
plt.savefig(newpath+r'\strain.jpg')

# FFT setup
from numpy.fft import fft, ifft
x = np.array(system.table.data['ll'])
fft_ll = fft(x)
N = len(fft_ll)
n = np.arange(N)
T = tfinal

tt = np.max(system.table.data['t'])
freq = n/T

# Plot the amplitude spectrum
plt.stem(freq, np.abs(fft_ll)/N, 'b', \
         markerfmt=" ", basefmt="-b")
plt.xlabel('f (Hz)')
plt.ylabel('Magnitude')
plt.xlim(0,5.5*source_freq)
plt.ylim(0,1)
plt.title('Model Amplitude Spectrum')
plt.savefig(newpath+r'\amp_spectrum.jpg')
