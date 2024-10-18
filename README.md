# Magnetostriction
A collection of code to describe magnetostrictive responses of ferromagnetic materials.

The goal of this repository is to determine what the strain signal will look like from a magnetostrictive material when positioned inside the sheathing of a fiber optic cable and placed in an alternating magnetic field. The files in this repository include:

<ul>
  <li>Time_Varying_Field_with_Magnetostriction.ipynb: Jupyter notebook to show the process of running a simulation using OOMMF.</li>
  <li>Time_Varying_Field_with_Magnetostriction.py: Runs one simulation with specified input parameters and outputs results to a specified folder.</li>
  <li>Time_Varying_Field_with_Magnetostriction_multi_sim.py: Runs a specified amount of simulations based on input parameters and incorporates temperature to determine outputs.</li>
  <li>PowerSpectrumPlot.ipynb: Jupyter notebook to plot the results from Time_Varying_Field_with_Magnetostriction_multi_sim.py. This file is not very friendly to use since the results of each simulation have different properties. I will document this file more.</li>
</ul>

Future additions to this repository will include descriptions of the magnetostrictive response of: 

<ul>
  <li>an anisotropic wire;</li> 
  <li>an anisotropic wire including eddy-current responses;</li>
  <li>a wire that is not in a uniform spatial magnetic field;</li>
  <li>and more...</li>
</ul>

If you are in this group, please feel free to contribute any code and/or literature that might provide insight into these magnetostrictive responses.
