%% Non-linear Analysis
% This is a look into the non-linearity of the system and what it means in
% terms of chaos theory. Hopefully this will provide some insight into the
% system's behavior.
clear;
close all

%% Plots to generate
% state space, delta theta, and poincare
state_plots = 0;

% Amplitude plots for harmonics at different H0
amp_plots = 1; 
HH0 = 1; % kA/m, source amplitude range of interest

% H0 varying with theta0
initial_conditions = 0; 
theta_i = 0:pi/10:pi; % radians, range of initial theta values 

% Plot harmonic amp with different f
f_plots = 0;
ff = 1:1:100; % Hz, range of frequency
omega_vec = 2*pi*ff; % rad/s, angular frequency

%% Adjustable Parameters
% Material Properties
Ms = 4.908e2; % kA/m, Saturation magnetization of Ni at 298K
L0 = -34; % ppm, saturation magnetostriction of Ni
lambda = 4.5e3; % adjustable damping parameter

% Source Properties, H = H0sin(omega*t-phi)
f = 100; % Hz, source frequency: omega = 2pi*f
H0 = 1; % kA/M, source amplitude
phi = 0; % radians, source phase - if set to -pi/2, the source becomes cos
         
% Initial Conditions and FFT Parameters
theta0 = [2*pi/3 2*pi/3.01]; % radians, initial theta positions. These should
                         % be close together to examine plots.
t_0 = 1; % s, time where theta and strain plots begin
t_f = 1.1; % s, final time
fs = 1e6; % Hz, FFT sampling frequency

%% Constants and some Calculations
g = 2; % spectroscopic splitting factor for e- spin
e = 1.6e-19; % C, e- charge
c = 3e8; % m/s, speed of light
me = 9.1e-31; % kg, e- mass

gamma = g*e/(2*me*c); % constant related to angular momentum
alpha = lambda/(gamma*Ms); % damping constant
omega = 2*pi*f; % Radians/s source frequency

%% ODE Solver (Runge-Kutta Method)
tspan = 0:1/fs:t_f-1/fs; % s, time span of interest
H = H0*sin(omega*tspan-phi); % magnetic field signal

t = zeros(length(tspan),length(theta0));
theta = zeros(length(tspan),length(theta0));
if state_plots == 1

for i = 1:length(theta0)
[t(:,i),theta(:,i)] = ...
    ode45(@(t,theta) LLG_2D(t,theta,gamma,alpha,H0,omega,phi,1,0),tspan,theta0(i));
end

% Solve for dthetadt and d2thetadt2 for analysis
dthetadt = LLG_2D(t,theta,gamma,alpha,H0,omega,phi,1,0);
d2thetadt2 = ...
    -gamma*(1+alpha^2)/(1+alpha)*H0*(omega*cos(omega*t-phi).*sin(theta)...
    + dthetadt.*sin(omega*t-phi).*cos(theta));

%% State Space
% Plot to see the dynamics in phase space
theta_0_pi = wrapToPi(theta);
t_initial = find(t(:,1)>=t_0);

figure(1)
subplot(2,1,1)
plot(theta(t_initial,1),dthetadt(t_initial,1),'-')
xlabel('\theta')
ylabel('$\dot{\theta}$','Interpreter','latex')
legend(strcat('\theta_0=',string(theta0(1))))
subplot(2,1,2)
plot(theta(t_initial,2),dthetadt(t_initial,2),'-',color=[0.9290 0.6940 0.1250])
xlabel('\theta')
ylabel('$\dot{\theta}$','Interpreter','latex')
legend(strcat('\theta_0=',string(theta0(2))))
sgtitle('State Space')

%% Difference in Thetas
delta = theta(:,2) - theta(:,1);
logdelta = log10(abs(delta));
dll = 3/2*L0.*(cos(theta(:,2)).^2-1/3)-3/2*L0.*(cos(theta(:,1)).^2-1/3); % ppm, homogeneous strain response

figure(2)
subplot(3,1,1)
plot(t(:,1),delta)
ylabel('\Delta\theta','FontSize',12)
title('Difference in Magnetic Moment Direction')
subplot(3,1,2)
plot(t(:,1),logdelta)
ylabel('log|\Delta\theta|','FontSize',12)
subplot(3,1,3)
plot(t(:,1),dll)
ylabel('\Delta\lambda_{\theta}','FontSize',12)
xlabel('time (s)','FontSize',12)
title('Difference in Strain')

%% Angular Acceleration State Space

figure(3)
subplot(2,1,1)
plot(dthetadt(t_initial,1),d2thetadt2(t_initial,1),'-')
xlabel('$\dot{\theta}$','Interpreter','latex')
ylabel('$\ddot{\theta}$','Interpreter','latex')
legend(strcat('\theta_0=',string(theta0(1))))
subplot(2,1,2)
plot(dthetadt(t_initial,2),d2thetadt2(t_initial,2),'-',color=[0.9290 0.6940 0.1250])
xlabel('$\dot{\theta}$','Interpreter','latex')
ylabel('$\ddot{\theta}$','Interpreter','latex')
legend(strcat('\theta_0=',string(theta0(2))))
sgtitle('Angular Acceleration State Space')

%% FFT Difference
ll = 3/2*L0.*(cos(theta).^2-1/3); % ppm, homogeneous strain response
y = fft(ll);

N = fs*t_f; %  number of samples
fq = fs*(0:N/2-1)/N; % convert to frequency domain
y = abs(y(1:N/2,:)/(N/2));

figure(10)
plot(fq,y(:,1),'-',fq,y(:,2),'-.')
title('FFT')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend(strcat('\theta_0=',string(theta0(1))),strcat('\theta_0=',string(theta0(2))))
xlim([0 20*f])



end

%% Frequency Variation with Amplitude
if amp_plots == true

%tH = zeros(length(t),length(HH0));
%thetaH = zeros(length(t),length(HH0));
dthetadtH = zeros(length(tspan),length(HH0));
for i = 1:length(HH0)
[tH(:,i),thetaH(:,i)] = ...
    ode45(@(t,theta) LLG_2D(t,theta,gamma,alpha,HH0(i),omega,phi,1,0),tspan,theta0(1));
    dthetadtH(:,i) = LLG_2D(tH(:,i),thetaH(:,i),gamma,alpha,HH0(i),omega,phi,1,0);
end

ll = 3/2*L0.*(cos(thetaH).^2-1/3); % ppm, homogeneous strain response
for i =1:length(HH0)
y(:,i) = fft(ll(:,i));
end

N = fs*t_f; %  number of samples
fq = fs*(0:N/2-1)/N; % convert to frequency domain
x_even = f*[4 8 12 16 18 20];
x_odd = f*[3 7 11 15 17 19];
x_tot = f*[1 2 3 4 5 6 7 8 0 10 11 12];
[~,fq_even] = ismember(x_even,fq);
[~,fq_odd] = ismember(x_odd,fq);
[~,fq_tot] = ismember(x_tot,fq);
y = abs(y(1:N/2,:)/(N/2));

figure(4)
plot(HH0,y(fq_tot(1:4),:),'.')
title('Frequency Magnitudes with Increasing Source Amplitude')
xlabel('Source Amplitude (kA/m)')
ylabel('Magnitude')
lgd = legend(strcat(string(x_tot(1:4)),' Hz'));

% figure(5)
% subplot(2,1,1)
% plot(HH0,y(fq_odd,:),'.')
% legend(strcat(string(x_odd),' Hz'))
% title('Odd Harmonics')
% ylabel('Magnitude')
% subplot(2,1,2)
% plot(HH0,y(fq_even,:),'.')
% legend(strcat(string(x_even),' Hz'))
% ylabel('Magnitude')
% xlabel('Source Amplitude (kA/m)')
% title('Even Harmonics')

% Bifurcation Diagram
syms theta t
t_use = tH(tH >= 1.09);
theta_use = thetaH(tH >= 1.09);
sol_vec = [];
sol1 = 0;
% for i= 1:length(t_use)
%     f = -gamma*(1+alpha^2)/(1+alpha)*H0.*sin(omega.*t-phi).*sin(theta);
%     sol = vpasolve(f,[theta, t],[theta_use(i), t_use(i)]);
% 
%     if vpa(sol.theta) ~= sol1
%         sol_vec = [sol_vec, vpa(sol.theta)];
%     end
% 
%     sol1 = vpa(sol.theta);
% end



figure(15)
plot(H0,sol_vec ,'.')
ylabel('\theta')
xlabel('H_0')
end



%% Variation of Frequency with Initial Conditions
if initial_conditions == 1

for i = 1:length(theta_i)
[t_th(:,i),theta_th(:,i)] = ...
    ode45(@(t,theta) LLG_2D(t,theta,gamma,alpha,H0,omega,phi),tspan,theta_i(i));
end

ll = 3/2*L0.*(cos(theta_th).^2-1/3); % ppm, homogeneous strain response
for i =1:length(theta_i)
y(:,i) = fft(ll(:,i));
end

N = fs*t_f; %  number of samples
fq = fs*(0:N/2-1)/N; % convert to frequency domain
x_even = f*[4 8 12 16 18 20];
x_odd = f*[3 7 11 15 17 19];
x_tot = f*[1 2 3 4 5 6 7 8 0 10 11 12];
[~,fq_even] = ismember(x_even,fq);
[~,fq_odd] = ismember(x_odd,fq);
[~,fq_tot] = ismember(x_tot,fq);
y = abs(y(1:N/2,:)/(N/2));

figure(6)
plot(theta_i,y(fq_tot(1:4),:),'.')
title('Frequency Magnitudes with \theta_0')
xlabel('\theta_0')
ylabel('Magnitude')
lgd = legend(strcat(string(x_tot(1:4)),' Hz'));

figure(7)
subplot(2,1,1)
plot(theta_i,y(fq_odd,:),'.')
legend(strcat(string(x_odd),' Hz'))
title('Odd Harmonics')
xlabel('\theta_0')
ylabel('Magnitude')
subplot(2,1,2)
plot(theta_i,y(fq_even,:),'.')
legend(strcat(string(x_even),' Hz'))
xlabel('\theta_0')
ylabel('Magnitude')
title('Even Harmonics')

end

%% Response Magnitude Frequency Dependence
if f_plots == true

%t_f = zeros(length(t),length(HH0));
%thetaf = zeros(length(t),length(HH0));
for i = 1:length(ff)
[t_ff(:,i),thetaf(:,i)] = ...
    ode45(@(t,theta) LLG_2D(t,theta,gamma,alpha,H0,omega_vec(i),phi),tspan,theta0(1));
end

ll = 3/2*L0.*(cos(thetaf).^2-1/3); % ppm, homogeneous strain response
for i =1:length(ff)
y(:,i) = fft(ll(:,i));
end

N = fs*t_f; %  number of samples
fq = fs*(0:N/2-1)/N; % convert to frequency domain

for j = 1:length(ff)
    x_even(j,:) = ff(j)*[4 8 12 16 18 20];
    x_odd(j,:) = ff(j)*[3 7 11 15 17 19];
    x_tot(j,:) = ff(j)*[1 2 3 4 5 6 7 8 0 10 11 12];
end
[~,fq_even] = ismember(x_even,fq);
[~,fq_odd] = ismember(x_odd,fq);
[~,fq_tot] = ismember(x_tot,fq);
y = abs(y(1:N/2,:)/(N/2));

figure(8)
plot(ff,y(fq_tot(1:2),:),'.')
title('Response Amplitudes with Source Frequency')
xlabel('Source Frequency (Hz)')
ylabel('Response Magnitude')
lgd = legend('Driving f','Double f');

% figure(9)
% subplot(2,1,1)
% plot(ff,y(fq_odd,:),'.')
% legend(strcat(string(x_odd),' Hz'))
% title('Odd Harmonics')
% ylabel('Magnitude')
% subplot(2,1,2)
% plot(ff,y(fq_even,:),'.')
% legend(strcat(string(x_even),' Hz'))
% ylabel('Magnitude')
% xlabel('Source Amplitude (kA/m)')
% title('Even Harmonics')

end
