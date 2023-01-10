function [tt, zeta] = ODESolver_asymp(deltat,T,n_ed, z, V, h,L)
% this file is used to solve the ODEs with symmetric salt
% if the system is asymmetric, the differential capacitance needs to be
% solved at each time step.......complicated in this case
tspan = linspace(0,T,floor(T/deltat)+1);
y0 = zeros(1,n_ed);
opts = odeset('Mass',@(t,q) mass_asymp(t,q,z,n_ed));
[tt,zeta] = ode45(@(t,q) fODE_asymp(t,q,z,n_ed,V,h,L),tspan,y0,opts);