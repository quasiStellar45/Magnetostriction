#!/usr/bin/env python
# coding: utf-8

# Import packages
import oommfc as mc
import discretisedfield as df
import micromagneticmodel as mm
import numpy as np
import matplotlib.pyplot as plt
import math
import random
import os
import pandas as pd
import scipy as sp

# define the number label of the run here, this will create a new folder for the outputs
newpath = r'C:\Users\Tomas\Desktop\Martin Group\3D_results\CWPpresentation_' # define model location and name here
# Run Number
num = 1        # initial number of run to start at
num_runs = 1  # number of simulations to run
# Time setup
tfinal = .5          # Final time of simulation
tstep =1e-5        # timestep in seconds
ntot = int(tfinal/tstep)         # total number of time steps
# Material properties
A_ex = 15e-12       # exchange constant
L1 = -24            # Strain constant 1 (<100> cubic direction) in ppm
L2 = -48            # Strain constant 2 (<111> cubic direction) in ppm
K1 = -0.5e5/1.257e2       # anisotropy constant 1 J/m^3: K1 < 0 indicates hard axis, K1 > 0 easy axis
K2 = -0.2e5/1.257e2       # anisotropy constant 2 J/m^3, not callable in this code 
u_11 = (1, 0, 0)    # K1 axis one: will be hard axis if K1 < 0, easy axis if K1 > 0
u_12 = (0, 1, 0)    # K1 axis 2: just needs to be orthogonal to u_11
damping = .75         # damping constant (alpha)
Ms = 5.35e5          # magnetization saturation (A/m)
Tc = 631            # K, Curie temperature
M0 = 5.2e5          # Ms at 0K (A/m)
# Alternating source properties
source_freq = 100   # source frequency (Hz)
H_app = (0, 0, 1.3e3)# applied field (A/m) in  (x,y,z) coordinates
source_type = 'sin' # waveform of source
# Strain measurement axis
B = np.array([0,0,1]) # ([x,y,z])
# Geometry of object
r = 1e-9            # material radius, m (this is only a radius if making a cylinder)
L = 2*r            # material length, m
d = 1e-9            # cell dimensions, m
# Temperature modeling
Temp = 294               # Temperature in K
'''
t = Temp/Tc              # fractional temperature

# function to determine magnetization saturation
def fn1(m):
    return m - np.tanh(m/t)

m = sp.optimize.brentq(fn1,0.001,2) # fractional mag sat 
Ms = M0*m                            # saturation at given temperature

# Sigmoid function definition
def sigmoid(x):
    return 0.12677279-0.10035642 / (1 + np.exp(-0.0851499*(x-69.35616788)))

# Apply the sigmoid function to temperature
damping = sigmoid(Temp)
'''

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
    
def cube_fun(pos):
    """Function to set values for a 4x4 cube

    """
    i_states = {'state_1':(1,0,1),'state_2':(-1,0,-1),'state_3':(0,1,0),'state_4':(0,-1,0)}
    x, y, z = pos
    if x > 0 and y > 0 and z > 0: 
        return i_states['state_1']
    if x > 0 and y < 0 and z > 0: 
        return i_states['state_3']
    if x < 0 and y < 0 and z > 0:
        return i_states['state_2']
    if x < 0 and y > 0 and z > 0:
        return i_states['state_4']
    if x > 0 and y > 0 and z < 0:
        return i_states['state_3']
    if x > 0 and y < 0 and z < 0:
        return i_states['state_2']
    if x < 0 and y < 0 and z < 0:
        return i_states['state_4']
    if x < 0 and y > 0 and z < 0:
        return i_states['state_1']
    
for j in range(0,num_runs):
    #display run mumber
    print("Run Number: ",str(num))
    newpath = newpath + str(num) # appends the number (num) of the run to the file name
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
    L2 = {str(L2)}                     # Strain constant 2 in ppm
    K1 = {str(K1)}                     # anisotropy constant 1
    K2 = {str(K2)}                     # anisotropy constant 2
    u_easy = {str(u_11)}               # axis 1 for K1
    u_hard = {str(u_12)}               # axis 2 for K1
    damping = {str(damping)}           # damping constant (alpha)
    Ms = {str(Ms)}                     # magnetization saturation (A/m)
    Tc = {str(Tc)}            # K, Curie temperature
    M0 = {str(M0)}          # Ms at 0K (A/m)
    # Alternating source properties
    source_freq = {str(source_freq)}   # source frequency (Hz)
    H_app = {str(H_app)}               # applied field (A/m)
    source_type = {source_type}        # waveform of source
    # Strain measurement axis
    B = {str(B)}
    # Geometry of object
    r = {str(r)}                       # material radius (m)
    L = {str(L)}                       # material length (m)
    d = {str(d)}                       # cell dimensions (m)
    # Temperature modeling
    T = {str(Temp)}               # Temperature in K
    """
    )

    f.close()

    # Define initial mesh
    p1 = (-r, -r, -L/2)                        # Starting point
    p2 = (r, r, L/2)                           # Ending point
    cell = (d, d, d)                           # Cell size
    region = df.Region(p1=p1, p2=p2)           # Define the region
    mesh = df.Mesh(region=region, cell=cell)   # Create the mesh

    # Initial mesh
    mesh.k3d()

    # Define the system name
    system = mm.System(name='time_dependent_field')

    # Define system energy
    system.energy = (mm.Zeeman(H=H_app, func=source_type, f=source_freq, t0=0) 
                    + mm.Demag()
                    + mm.Exchange(A=A_ex)
                    + mm.CubicAnisotropy(K=K1, u1=u_11, u2=u_12))

    # Define system dynamics
    system.dynamics = mm.Precession(gamma0=mm.consts.gamma0) + mm.Damping(alpha=damping)

    # Define the field element
    system.m = df.Field(mesh, dim=3, value=cube_fun, norm=Ms)

    # Plot the initial magnetization of one cross section that is one cell thick
    system.m.plane('z').mpl()
    plt.savefig(newpath+r'\initial_mag.jpg')

    ev = mc.RungeKuttaEvolver(min_timestep=tstep) # Define the evolver as Runge Kutta method
    td = mc.TimeDriver(evolver=ev)                # Setup time driver with evolver input
    td.drive(system, n=ntot, t=tfinal, verbose=2) # Drive the system

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
    ll_rate = [0]
    # calculate strain rate by numerical differentiation over one time step
    for i in range(0,len(ll)-1):
        ll_rate_i = (ll[i+1]-ll[i])/(system.table.data['t'][1]-system.table.data['t'][0])
        ll_rate.append(ll_rate_i)

    # Define a new column in the data for the strain
    system.table.data['ll'] = ll
    system.table.data['ll_rate'] = ll_rate

    # Plot the strain response
    system.table.mpl(y=['ll'])
    plt.legend().remove()
    plt.title('Model Strain Response')
    plt.ylabel(r'$\lambda$ (ppm)')
    plt.savefig(newpath+r'\strain.jpg')

    # Plot the strain response zoomed
    #system.table.mpl(y=['ll'])
    #plt.legend().remove()
    #plt.title('Model Strain Response')
    #plt.ylabel(r'$\lambda$ (ppm)')
    #plt.xlim([28,30])
    #plt.ylim([7.5875,7.595])
    #plt.savefig(newpath+r'\strain_zoom.jpg')

    # plot strain rate response
    system.table.mpl(y=['ll_rate'])
    plt.legend().remove()
    plt.title('Model Strain Rate')
    plt.ylabel(r'$\frac{d \lambda}{dt}$ (ppm/s)')
    plt.savefig(newpath+r'\strain_rate.jpg')

    # zoomed strain rate response
    #system.table.mpl(y=['ll_rate'])
    #plt.legend().remove()
    #plt.title('Model Strain Rate')
    #plt.ylabel(r'$\frac{d \lambda}{dt}$ (ppm/s)')
    #plt.ylim([-500,500])
    #plt.xlim([29.8,30])
    #plt.savefig(newpath+r'\strain_rate_zoom.jpg')

    # FFT setup
    from scipy.fft import fft, fftfreq
    T = system.table.data['t'][1]-system.table.data['t'][0]
    y = np.array(system.table.data['ll'])

    fft_ll = fft(y)
    N = len(y)
    x = np.linspace(0.0, N*T, N, endpoint=False)

    xf = fftfreq(N, T)[:N//2]

    PS = 2*np.abs(fft_ll[0:N//2])/N

    # Plot the amplitude spectrum
    #plt.plot(xf, PS)
    #plt.xlabel('f (Hz)')
    #plt.ylabel('Magnitude')
    #plt.xlim(0,5*source_freq)
    #plt.ylim(0,0.001)
    #plt.grid()
    #plt.title('Model Amplitude Spectrum')
    #plt.savefig(newpath+r'\amp_spectrum.jpg')

    # Export Data
    system.table.data.to_csv(newpath+r'/data.csv')
    dff = pd.DataFrame({'frequency':xf, 'Amplitude':PS})
    dff.to_csv(newpath+r'/Power_Spectrum.csv')

    # FFT setup
    T = system.table.data['t'][1]-system.table.data['t'][0]
    y = np.array(system.table.data['ll'][system.table.data['t']>.1])

    fft_ll = fft(y)
    N = len(fft_ll)
    x = np.linspace(0.0, N*T, N, endpoint=False)

    xf = fftfreq(N, T)[:N//2]

    PS = 2*np.abs(fft_ll[0:N//2])/N

    # Record the amplitude at 100Hz, 200Hz, and harmonics up to 1000Hz
    peak_f = [100,200,300,400,500,600,700,800,900,1000]
    peaks = [PS[40],PS[80],PS[120],PS[160],PS[200],PS[240],PS[280],PS[320],PS[360],PS[400]]
    df_peak = pd.DataFrame({'frequency':peak_f, 'Amplitude':peaks})
    df_peak.to_csv(newpath+r'/peak_values.csv')
    # Plot the amplitude spectrum
    #plt.plot(xf, PS)
    #plt.grid()
    #plt.xlabel('f (Hz)')
    #plt.ylabel('Magnitude')
    #plt.xlim(0,5*source_freq)

    #plt.ylim(0,.2)
    #plt.title('Filtered Model Amplitude Spectrum')
    #plt.savefig(newpath+r'\amp_spectrum_filtered.jpg')

    # FFT setup for strain rate
    T = system.table.data['t'][1]-system.table.data['t'][0]
    y = np.array(system.table.data['ll_rate'][system.table.data['t']>.1])

    fft_ll = fft(y)
    N = len(fft_ll)
    x = np.linspace(0.0, N*T, N, endpoint=False)

    xf_rate = fftfreq(N, T)[:N//2]

    PS_rate = 2*np.abs(fft_ll[0:N//2])/N

    # Plot the amplitude spectrum of the strain rate
    #plt.plot(xf_rate, PS_rate)
    #plt.xlabel('f (Hz)')
    #plt.ylabel('Magnitude')
    #plt.xlim(0,5*source_freq)
    #plt.ylim(0,0.01)
    #plt.grid()
    #plt.title('Model Strain Rate Amplitude Spectrum')
    #plt.savefig(newpath+r'\amp_spectrum_filtered_strain_rate.jpg')

    df_rate = pd.DataFrame({'frequency':xf_rate, 'Amplitude':PS_rate})
    df_rate.to_csv(newpath+r'/Rate_Power_Spectrum.csv')

    # edit damping constant for next run
    #damping += 1         # damping constant (alpha)
    # edit magnetic field strength for next run
    # H_app = tuple(np.array(H_app) + np.array((0,0,.01)))
    # edit temperature for next run
    #Temp += 1
    #edit run number for next run
    num += 1