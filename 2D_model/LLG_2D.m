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
