%% Model of Magnetic Moment Motion
clear;
close all

%% Adjustable Parameters
f = 100; % Hz, source frequency
L0 = -34; % ppm, saturation magnetostriction of Ni
H0 = 1;% A/M, source amplitude
Ms = 0.124; % Saturation magnetization of Ni
lambda = 1e-6; % adjustable damping parameter
theta0 = pi/3; % radians, initial theta position
t_f = 120; % s, final time
fs = 10000; % Hz, FFT sampling frequency

%% Run Code
g = 2; % spectroscopic splitting factor for e- spin
e = 1.6e-19; % C, e- charge
c = 3e8; % m/s, speed of light
me = 9.1e-31; % kg, e- mass

gamma = g*e/(2*me*c); % constant related to angular momentum
alpha = lambda/(gamma*Ms); % damping constant
omega = 2*pi*f; % Radians/s source frequency
tspan = 0:1/fs:t_f-1/fs; % s, time span of interest

[t,theta] = ode45(@(t,theta) fun(t,theta,gamma,alpha,H0,omega),tspan,theta0);


figure(1)
plot(t,theta,'-')
xlabel('Time (s)');
ylabel('\theta (radians)');
title('Magnetic Moment Angle')

% Strain
ll = 3/2*L0.*(cos(theta).^2-1/3); % ppm, homogeneous strain response

figure(2)
plot(t,ll,'-')
xlabel('Time (s)');
ylabel('\lambda (ppm)');
title('Strain Response')

% Fourier Transform
y = fft(ll);
N = fs*t_f; %  number of samples
f = fs*(0:N/2-1)/N; % convert to frequency domain

figure(3)
plot(f,abs(y(1:N/2)/(N/2)))
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Fourier Transform')
ylim([0 5])

function dthetadt = fun(t,theta,gamma,alpha,H0,omega)
dthetadt = gamma*(1+alpha)/(1+alpha^2)*H0.*sin(omega.*t).*sin(theta);
end