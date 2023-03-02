# Magnetostriction
A collection of code to describe magnetostrictive responses of ferromagnetic materials.

The goal of this repository is to determine what the strain signal will look like from a magnetostrictive material when positioned inside the sheathing of a fiber optic cable and placed in an alternating magnetic field. The files in this repository include:

<ul>
  <li>Domain_Wall.m : This file attempts to describe the motion of domain walls within a ferromagnetic material when placed in an alternating magnetic field. The formula used for the domain wall motion is $m\frac{d^2x}{dt^2}+\beta \frac{dx}{dt}+\alpha x = 2M_sH$.</li>
  <li>Isotropic_2D_Magnetostrictive_Response.m : This file uses the Landau-Lifshitz-Gilbert equation, <br> $\frac{d\vec{M}}{dt}=\gamma(\vec{M}\times \vec{H})-\frac{\alpha}{M}\left(\vec{M}\times\frac{d\vec{M}}{dt}\right)+\gamma \alpha^2(\vec{M}\times \vec{H})$, <br> to describe the motion of the magnetic dipole of a single domain. The crystal lattice is assumed to experience an isotropic magnetostrictive response, thus simplifying the magnetostrictive formula to $\lambda_{\theta}=\frac{3}{2}\lambda_s(cos^2\theta - \frac{1}{3})$. The file "2D Isotropic Magnetostriction Model.pptx" has some more information on the output of this code. </li>
</ul>

Future additions to this repository will include descriptions of the magnetostrictive response of: 

<ul>
  <li>multiple domains; </li>
  <li>a 3D system; </li>
  <li>a anisotropic domain;</li> 
  <li>an anisotropic wire;</li> 
  <li>an anisotropic wire including eddy-current responses;</li>
  <li>a wire that is not in a uniform spatial magnetic field;</li>
  <li>and more...</li>
</ul>

If you are in this group, please feel free to contribute any code and/or literature that might provide insight into these magnetostrictive responses.
