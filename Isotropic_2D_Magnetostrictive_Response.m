%% Model of Isotropic, 2D, Single Domain Magnetostricitive Response 
% The assumptions are that the material is homogeneous and isotropic, 
% there are 2 dimensions, there is one domain, the magnetic field is 
% uniform in space, the strain measurement direction is in the same 
% direction as the magnetic field, and the magnetic moment torque response 
% is instantaneous. The source is sinusoidal with changeable phase and 
% frequency. The damping parameter is unknown and can be altered to 
% examine changes in the response.
clear;
close all

%% Which Plots to Generate - enter true or false (1 or 0)
dynamics_plots = 1;   % theta dynamics plots
strain_plot = 0;      % strain plot
FFT_plot = 1;         % FFT plot
hysteresis_plots = 0; % hysteresis plots
run_large = 0;        % large scale animation
run_small = 0;        % small scale animation
                      % - both animations cannot run at the same time
                      % - make sure to check parameters in animations
                      %   section before running
%% Adjustable Parameters
% Material Properties
Ms = 4.908e2; % kA/m, Saturation magnetization of Ni at 298K
L0 = -34;     % ppm, saturation magnetostriction of Ni
lambda = 1;   % adjustable damping parameter

% Source Properties, H = H0sin(omega*t-phi)+H_bias or
% H0square(omega*t)+Hbias
HSource = 1;  % set to 1 for a sin source, 0 for a square wave. To 
              % run the square wave you will need MATLAB's signal
              % processing package.
f = 100;      % Hz, source frequency: omega = 2pi*f 
H0 = 8;      % kA/m, source amplitude
phi = 0;      % radians, source phase; if set to -pi/2, the source becomes 
              % cosine.
H_bias = 0;   % Add a bias to the source signal. If =0, no bias is included

% Initial Conditions and FFT Parameters
theta0 = pi/3; % radians, initial theta position
t_0 = 0;       % s, time where theta and strain plots begin
t_f = 2;       % s, final time
fs = 1e5;      % Hz, FFT sampling frequency
sety = 10;     % yaxis upper limit of the FFT

%% Animations - check these parameters before running the animations
% Large Scale Animation - displays all motion of the magnetic moment
step_large = 10;     % iteration step of animation
Hfactor_large = 5;   % factor to increase the magnetic field vector by. 
                       % This value should be H_factor*H0 <= Ms.
t_frac_large = 100;   % What fraction of the total time do you want to run?

% Small Scale Animation - displays small scale motion of magnetic moment
% when H is in the theta = 0 direction. Run the large scale animation first
% to see if there are small oscillations. If there are, run this animation
% to increase the size of those small scale oscillations.
step_small = 1;       % iteration step of animation
Hfactor_small = 5;    % factor to increase the magnetic field vector by
t_frac_small = 200;    % What fraction of the total time do you want to run?
theta_limit_0 = 5e-10; % upper limit of thetas of interest when the magnetic 
                      % field is pointing in the direction of theta = 0. 
                      % Note that this may need to be on the order of 
                      % 10^-10 for high H0 (>10 kA/m).
theta_factor_0 = 1e5; % factor to increase theta by when the magnetic 
                      % field is pointing in the theta = 0 direction. 
                      % Note that this may need to be on the order of 10^5 
                      % for high H0 and low theta_limit_0.
                    
%% Constants and some Calculations
g = 2; % spectroscopic splitting factor for e- spin
e = 1.6e-19;  % C, e- charge
c = 3e8;      % m/s, speed of light
me = 9.1e-31; % kg, e- mass
k = pi*4e-7;  % H/m

gamma = g*e/(2*me*c);      % constant related to angular momentum
alpha = lambda/(gamma*Ms); % damping constant
omega = 2*pi*f;            % Radians/s source frequency

%% ODE Solver (Runge-Kutta Method)
tspan = 0:1/fs:t_f-1/fs;     % s, time span of interest
if HSource == 1
H = H0*sin(omega*tspan-phi)+H_bias; % magnetic field signal
else
H = H0*square(omega*tspan);
end

[t,theta] = ...
    ode45(@(t,theta) LLG_2D(t,theta,gamma,alpha,H0,omega,phi,HSource,H_bias),tspan,theta0);

%% Plot of theta
if HSource == 1
    j = 1;
else
    j = 0;
end

if dynamics_plots == true
figure(1)
subplot(2+j,1,1)
plot(t,theta,'-')
xlabel('Time (s)');
ylabel('\theta (radians)');
title('Magnetic Moment Angle')
xlim([t_0 t_f])
yyaxis right
plot(t,H,'r')
ylim([-10*H0 10*H0])
ylabel('H (kA/m)')
legend('\theta','H')

%% Plot of Angular Velocity
dthetadt = LLG_2D(t,theta,gamma,alpha,H0,omega,phi,HSource,H_bias); % Calc dthetadt 

subplot(2+j,1,2)
plot(t,dthetadt,'-')
xlabel('Time (s)');
ylabel('$\dot{\theta} (radians/s)$','Interpreter','latex');
title('Magnetic Moment Angular Velocity')
xlim([t_0 t_f])
yyaxis right
plot(t,H,'r')
ylim([-10*H0 10*H0])
ylabel('H (kA/m)')
vl = legend('$\dot{\theta}$','H');
set(vl, 'Interpreter', 'latex');

%% Plot of Angular Acceleration
% Calculate angular acceleration (only for sine wave)
if HSource == 1
d2thetadt2 = ...
    -gamma*(1+alpha^2)/(1+alpha)*H0*(omega*cos(omega*t-phi).*sin(theta)...
    + dthetadt.*(sin(omega*t-phi)+H_bias).*cos(theta));

subplot(3,1,3)
plot(t,d2thetadt2,'-')
xlabel('Time (s)');
ylabel('$\ddot{\theta} (radians/s^2)$','Interpreter','latex');
title('Magnetic Moment Angular Acceleration')
xlim([t_0 t_f])
yyaxis right
plot(t,H,'r')
ylim([-10*H0 10*H0])
ylabel('H (kA/m)')
al = legend('$\ddot{\theta}$','H');
set(al, 'Interpreter', 'latex');
end

end

%% Plot of Strain
ll = 3/2*L0.*(cos(theta).^2-1/3); % ppm, homogeneous strain response
if strain_plot == true

figure(2)
plot(t,ll,'-')
xlabel('Time (s)');
ylabel('\Lambda_{\theta} (ppm)');
title('Strain Response')
xlim([t_0 t_f])
yyaxis right
plot(t,H,'r')
ylim([-10*H0 10*H0])
ylabel('H (kA/m)')
legend('\Lambda_{\theta}','H')
end

%% Fourier Transform
if FFT_plot == true
y = fft(ll);
N = fs*t_f; %  number of samples
fq = fs*(0:N/2-1)/N; % convert to frequency domain

figure(3)
plot(fq,abs(y(1:N/2)/(N/2)))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Fourier Transform')
ylim([0 sety])
xlim([0 50*f])
end

%% Magnetization Hysteresis Plot
if hysteresis_plots == true
Mx = Ms*cos(theta);
My = Ms*sin(theta);

figure(4)
subplot(1,2,1)
plot(H,Mx)
xlabel('H (kA/m)')
ylabel('M_x (kA/m)')
grid on
subplot(1,2,2)
plot(H,My)
xlabel('H (kA/m)')
ylabel('M_y (kA/m)')
grid on

sgtitle('Magnetization Hysteresis Plots')

%% Strain Hysteresis Plot
figure(5)
plot(H,ll)
xlabel('H (kA/m)')
ylabel('\Lambda_{\theta} (ppm)')
title('Strain Hysteresis Plot')
end

%% Motion Animation
% The magnetic moment is blue and the magnetic field is red. Hfactor
% mutliplies the magnetic field by a specified constant.
if run_large == true
    for i=1:step_large:length(t)/t_frac_large
        figure(6)
        polarplot(theta(i)*(0:1),0:1,'b-',[0 0],Hfactor_large*[0 H(i)/Ms],'r-')
        title(strcat('Time:',{' '},string(t(i)),'s'))
    end
end

%% Small Motion Animation
% This animation displays the small angle motion on a polar plot and
% amplifies it so it can be seen better.
if run_small == true
    for i = 1:step_small:length(t)/t_frac_small
        if (abs(1-cos(theta(i))) < theta_limit_0)
            figure(7)
            polarplot(theta_factor_0*theta(i)*(0:1),0:1,'b-',[0 0],Hfactor_small*[0 H(i)/Ms],'r-')
            title(strcat('Time:',{' '},string(t(i)),'s'))
        else
            figure(7)
            polarplot(0,1,[0 0],Hfactor_small*[0 H(i)/Ms],'r-')
            title(strcat('Time:',{' '},string(t(i)),'s'))
        end
    end
end

%% Equation of Motion
function dthetadt = LLG_2D(t,theta,gamma,alpha,H0,omega,phi,HSource,H_bias)
if HSource == 1
dthetadt = ...
    -gamma*(1+alpha^2)/(1+alpha)*(H0.*sin(omega.*t-phi)+H_bias).*sin(theta);
else
dthetadt = ...
    -gamma*(1+alpha^2)/(1+alpha)*(H0.*square(omega*t)+H_bias).*sin(theta);
end
end