%... this file is used to solve the eigenvalues
clear
clc
tic
%% ode system of size 2n*2n
%%
H = 0.5;% assume
L = 0.5;
N =6;%2n
eps = (H+L)/100;
x_ed = zeros(1, N);
n_0 =floor(N/2);
h = H/(n_0-1);
deltax = 1/2000;
V =[-2,2];
z = [1,-1];
T =20;70;100;
deltat= deltax*0.2;
alpha = z(1,2)^2*z(1,1) - z(1,1)^2*z(1,2);
% eps = h/25;;
b =zeros(1, N);
c = zeros(1, N-1);
b(1,1:N) =[ 1/h,2/h*ones(1,N/2-2),1/L/2+1/h,1/L/2+1/h,2/h*ones(1,N/2-2),1/h];
c(1,1:N-1) = [-1/h,-1/h*ones(1,N/2-2),-1/2/L,-1/h*ones(1,N/2-1)];%c(N)=0 not used
a = [-1/h*ones(1,N/2-1),-1/2/L,-1/h*ones(1,N/2-2),-1/h];
M =  sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c,a],N,N);
M0 = M*h;
C_star = C_file(N,deltat,T, z, V, h,L);
M_star = diag(C_star)\M0;
eig_value2n = sort(eig(full(M_star)))
tau = h/(alpha*eig_value2n(2,1))




