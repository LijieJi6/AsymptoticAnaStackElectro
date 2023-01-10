function [tt, zeta] = ODESolver2(deltat,T,len, z, V, h,L,T_real,y_initial)
% this file is used to solve the ODEs with symmetric salt
% if the system is asymmetric, the differential capacitance needs to be
% solved at each time step.......complicated in this case

tspan = linspace(T_real+deltat,T,((T-T_real)/deltat)+1);%T,floor(T/deltat));%;linspace(deltat,T,floor(T/deltat));
y0 = y_initial;%zeros(1,2*(floor(length(len)/2)+1));
opts = odeset('Mass',@(t,q) mass2(t,q,z,len));
[tt,zeta] = ode45(@(t,q) fODE(t,q,z,len,V,h,L),tspan,y0,opts);