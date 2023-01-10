function C_star = C_file(N,deltat,T, z, V, h,L)

[~, zeta_ode] = ODESolver_asymp(deltat,T, N,z, V, h,L);
save('zeta_ode.mat','zeta_ode')
q = zeta_ode(end,:);
C_star = sqrt(2)*power(abs(z(1,1)),3/2)*cosh(z(1,1)*q/2);
C_star(2:end-1) = 2*C_star(2:end-1);