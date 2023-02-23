clear; close all;

lam100 = -42e-6; % Nickel easy axis strain
Ms = 55.1; % magnetic saturation of Nickel, J T^-1 kg^-1
k = 4*pi*10^(-7); % H/m
m = 1e-9; % mass per unit area of domain wall
beta = 10e-5; % viscous damping parameter
delta = 100000; 
alpha = 6*delta*lam100*k; % force due to crystal imperfections

f = 1; % Hz, frequency of source
omega = 2*pi*f; % angular frequency of source
H0 = 1e-6; % T, source amplitude

[t,x] = ode45(@(t,x) domainmotion(t,x,Ms,H0,omega,beta,alpha,m),[0,20],[0;0]);

figure(1)
plot(t,x(:,1),'-')
xlabel('Time t');
ylabel('Solution x');

y = ones(size(x(:,1)));

function dxdt = domainmotion(t,x,Ms,H0,omega,beta,alpha,m)
    dxdt = [x(2); (2*Ms*H0*sin(omega*t)-beta*x(2)-alpha*x(1))/m];
end


