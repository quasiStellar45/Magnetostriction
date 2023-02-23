# Magnetostriction
A collection of code to describe magnetostrictive responses of ferromagnetic materials.

The goal of this repository is to determine what the strain signal will look like from a magnetostrictive material when positioned inside the sheathing of a fiber optic cable and placed in an alternating magnetic field. The files in this repository include:

\begin{itemize}
\item Domain_Wall.m : This file attempts to describe the motion of domain walls within a ferromagnetic material when placed in an alternating magnetic field. The formula used for the domain wall motion is $m\frac{d^2x}{dt^2}+\Beta \frac{dx}{dt}+\alpha x = 2M_sH$.
\item Moment_Motion : This file uses the Landau-Lifshitz-Gilbert equation to describe the motion of the magnetic dipole of a single domain. The crystal lattice is assumed to experience an isotropic magnetostrictive response, thus simplifying the magnetostrictive formula to $\lambda_{\theta}=\frac{3}{2}\lambda_s(cos^2\theta - \frac{1}{3})$.
\end{itemize}

Future additions to this repository will include descriptions of the magnetostrictive response of: 

\begin{itemize}
\item a anisotropic domain; 
\item an anisotropic wire; 
\item an anisotropic wire including eddy-current responses;
\item a wire that is not in a uniform spatial magnetic field;
\item and more...
\end{itemize}
