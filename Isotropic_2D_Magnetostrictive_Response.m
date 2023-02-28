%% Model of Magnetostricitive Response 
% The assumptions are that the material is isotropic, there are 2
% dimensions, there is one domain, the magnetic field is uniform in space,
% and the magnetic moment torque response is instantaneous. The source is
% sinusoidal with changeable phase and frequency. The damping parameter is
% unknown and can be altered to examine changes in the response.
clear;
close all

%% Adjustable Parameters
% Material Properties
Ms = 4.908e5; % A/m, Saturation magnetization of Ni at 298K
L0 = -34; % ppm, saturation magnetostriction of Ni
lambda = 1e-6; % adjustable damping parameter

% Source Properties
f = 1; % Hz, source frequency
H0 = 1e4; % A/M, source amplitude
phi = 0; % radians, source phase - if set to pi/2, the source becomes cos

% Initial Conditions and FFT Parameters
theta0 = pi/4; % radians, initial theta position
t_0 = 0; % s, time where theta and strain plots begin
t_f = 120; % s, final time
fs = 1e5; % Hz, FFT sampling frequency

%% Differential Equation Solver
g = 2; % spectroscopic splitting factor for e- spin
e = 1.6e-19; % C, e- charge
c = 3e8; % m/s, speed of light
me = 9.1e-31; % kg, e- mass

gamma = g*e/(2*me*c); % constant related to angular momentum
alpha = lambda/(gamma*Ms); % damping constant
omega = 2*pi*f; % Radians/s source frequency
tspan = 0:1/fs:t_f-1/fs; % s, time span of interest

H = H0*sin(omega*tspan+phi); % magnetic field signal

[t,theta] = ode45(@(t,theta) fun(t,theta,gamma,alpha,H0,omega,phi),tspan,theta0);

%% Plot of theta
figure(1)
plot(t,theta,'-')
xlabel('Time (s)');
ylabel('\theta (radians)');
title('Magnetic Moment Angle')
xlim([t_0 t_f])
yyaxis right
plot(t,H,'r')
ylim([-10*H0 10*H0])
ylabel('H (A/m)')
legend('\theta','H')

%% Converted plot of theta to radians from 0 to 2pi
conv_sin = sin(theta); % taking sin of theta for conversion 
conv_cos = cos(theta); % taking cosin of theta for conversion

% Converting the angles to values from -pi to pi
theta_rad = zeros(size(theta));
for i = 1:length(t)
    if (conv_sin(i) >= 0) && (conv_cos(i) > 0)
        theta_rad(i) = acos(conv_cos(i));
    elseif (conv_sin(i) > 0) && (conv_cos(i) <= 0)
        theta_rad(i) = asin(conv_sin(i))+pi/2;
    elseif (conv_sin(i) <= 0) && (conv_cos(i) < 0)
        theta_rad(i) = acos(conv_cos(i))-3*pi/2;
    elseif (conv_sin(i) < 0) && (conv_cos(i) >= 0)
        theta_rad(i) = asin(conv_sin(i));
    else
        print('error')
    end
end

figure(2)
plot(t,theta_rad)
xlabel('Time (s)');
ylabel('\theta (radians)');
title('Magnetic Moment Angle')
xlim([t_0 t_f])
yyaxis right
plot(t,H,'r')
hold on
plot(t,zeros(length(t),1))
ylim([-H0 10*H0])
ylabel('H (A/m)')
legend('\theta','H')
hold off

%% Plot of Strain
ll = 3/2*L0.*(cos(theta).^2-1/3); % ppm, homogeneous strain response

figure(3)
plot(t,ll,'-')
xlabel('Time (s)');
ylabel('\lambda_{\theta} (ppm)');
title('Strain Response')
xlim([t_0 t_f])
yyaxis right
plot(t,H,'r')
ylim([-10*H0 10*H0])
ylabel('H (A/m)')
legend('\theta','H')

%% Fourier Transform
y = fft(ll);
N = fs*t_f; %  number of samples
fq = fs*(0:N/2-1)/N; % convert to frequency domain

figure(4)
plot(fq,abs(y(1:N/2)/(N/2)))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Fourier Transform')
ylim([0 5])
xlim([0 100*f])


function dthetadt = fun(t,theta,gamma,alpha,H0,omega,phi)
dthetadt = gamma*(1+alpha^2)/(1+alpha)*H0.*sin(omega.*t+phi).*sin(theta);
end