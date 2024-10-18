# Two-dimensional model
This folder contains code used for two-dimensional modeling of magnetostriction. This code was created with MATLAB 2022b. The code files are:
<ul>
  <li>Domain_Wall.m : This file attempts to describe the motion of domain walls within a ferromagnetic material when placed in an alternating magnetic field. The formula used for the domain wall motion is $m\frac{d^2x}{dt^2}+\beta \frac{dx}{dt}+\alpha x = 2M_sH$.</li>
  <li>Isotropic_2D_Magnetostrictive_Response.m : This file uses the Landau-Lifshitz-Gilbert equation, <br> $\frac{d\vec{M}}{dt}=\gamma(\vec{M}\times \vec{H})-\frac{\alpha}{M}\left(\vec{M}\times\frac{d\vec{M}}{dt}\right)+\gamma \alpha^2(\vec{M}\times \vec{H})$, <br> to describe the motion of the magnetic dipole of a single domain. The crystal lattice is assumed to experience an isotropic magnetostrictive response, thus simplifying the magnetostrictive formula to $\lambda_{\theta}=\frac{3}{2}\lambda_s(cos^2\theta - \frac{1}{3})$. The file "2D Isotropic Magnetostriction Model.pptx" has some more information on the output of this code. </li>
  <li>Nonlinear_Analysis.m : This file contains code to perform nonlinear analysis on the two-dimensional model. Outputs include phase space diagrams, Poincare plots, and initial angle difference plots.</li>
  <li>LLG_2D.m : A function for numerically solving the two-dimensional LLG equation. This function is used in Isotropic_2D_Magnetostrictive_Response.m and Nonlinear_Analysis.m.</li>
</ul>
