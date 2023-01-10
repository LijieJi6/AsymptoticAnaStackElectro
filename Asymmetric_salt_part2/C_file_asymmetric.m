function [q,C_star] = C_file_asymmetric(N,deltat,T, z, V, h,L,T_real)

[~, zeta_ode1] = ODESolver1(deltat,N, z, V, h,L,T_real);
y_initial =zeta_ode1(end,:);
save('ode_zeta1.mat','zeta_ode1')
[~, zeta_ode2] = ODESolver_asymmetric(deltat,T, N,z, V, h,L,T_real,y_initial);
save('ode_zeta2.mat','zeta_ode2')
% save('zeta_ode.mat','zeta_ode')
q = zeta_ode2(end,:);
len0 = floor(N/2);
C11=(z(1,1)*z(1,2)*exp(-z(1,1)*q(1:len0)) -z(1,1)*z(1,2)*exp(-z(1,2)*q(1:len0)));
C12 = (-sqrt(2)*sqrt(-z(1,2)*exp(-z(1,1)*q(1:len0)) +z(1,1)*exp(-z(1,2)*q(1:len0))+z(1,2)-z(1,1)));
C21=(z(1,1)*z(1,2)*exp(-z(1,1)*q(len0+1:2*len0)) -z(1,1)*z(1,2)*exp(-z(1,2)*q(len0+1:2*len0)));
C22 = (sqrt(2)*sqrt(-z(1,2)*exp(-z(1,1)*q(len0+1:2*len0)) +z(1,1)*exp(-z(1,2)*q(len0+1:2*len0))+z(1,2)-z(1,1)));
C11(2:len0) =2*C11(2:len0);
C21(1:len0-1) = 2*C21(1:len0-1);
save('C.mat','C21','C11')
C_star = [C11./C12,C21./C22];