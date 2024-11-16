# Application models
This folder contains code used to model application areas of interest. To use this code, you will need to download the [SimPEG](https://simpeg.xyz/) github [repository](https://github.com/simpeg/simpeg). I ran these codes by creating my own folder in the github folder download called "my_tests". The code in this folder is:

<ul>
  <li>DMS_sim_seawater.py : This code is used to simulate seawater instrusion monitoring. The outputs include a TMI anomaly plot and result files.</li>
  <li>DMS_sim_tailings.py : This code is used to simulate mine drainage and lithium brine monitoring. The outputs include a TMI anomaly plot and result files.</li>
  <li>DMS_multi_contour_plot.py : This code is used to plot the contour along y=0m of the result files from DMS_sim_seawater.py and DMS_sim_tailings.py.</li>
  <li>DMS_inversion.py : This file is used to invert results from DMS_sim_seawater.py and DMS_sim_tailings.py. This code did not seem to work well for negative values of magnetic susceptibility.</li>
  <li>DMS_sim_seawater_borehole.py : This code is used to simulate seawater intrusion inside a borehole.</li>
  <li>DMS_sim_metals_borehole.py : This code is used to simulate metal ion boundaries inside a borehole.</li>
  <li>DMS_contour_plot_borehole.py : This code is used to plot the contour of TMI anomaly from DMS_sim_seawater_borehole.py.</li>
</ul>

