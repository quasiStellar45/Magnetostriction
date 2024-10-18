# Three-dimensional model
This folder contains code used for three-dimensional modeling of magnetostriction. To run this code, install [Ubermag](https://ubermag.github.io/installation.html), following the conda install instructions. The code in this folder is:

<ul>
  <li>Time_Varying_Field_with_Magnetostriction.ipynb: Jupyter notebook to show the process of running a simulation using OOMMF.</li>
  <li>Time_Varying_Field_with_Magnetostriction.py: Runs one simulation with specified input parameters and outputs results to a specified folder.</li>
  <li>Time_Varying_Field_with_Magnetostriction_multi_sim.py: Runs a specified amount of simulations based on input parameters and incorporates temperature to determine outputs.</li>
  <li>PowerSpectrumPlot.ipynb: Jupyter notebook to plot the results from Time_Varying_Field_with_Magnetostriction_multi_sim.py. This file is not very friendly to use since the results of each simulation have different properties. I will document this file more.</li>
</ul>
