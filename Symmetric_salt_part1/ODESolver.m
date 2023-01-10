function [tt, zeta] = ODESolver(deltat,T,len, z, V, h,L)
% this file is used to solve the ODEs with symmetric salt
% if the system is asymmetric, the differential capacitance needs to be
% solved at each time step.......complicated in this case

tspan = linspace(0,T,floor(T/deltat)+1);
y0 = zeros(1,floor(length(len)/2)+1);
opts = odeset('Mass',@(t,q) mass(t,q,z,len));
[tt,zeta] = ode45(@(t,q) fODE(t,q,z,len,V,h,L),tspan,y0,opts);