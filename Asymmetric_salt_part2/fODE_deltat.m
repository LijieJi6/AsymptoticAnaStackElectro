
function dydt = fODE_deltat(t, q, z, n_ed, V, h,L)
N = n_ed;% number of total electrodes 2n
alpha = z(1,2)^2*z(1,1) - z(1,1)^2*z(1,2);
if N==1
    M0 = 1;
    dydt = z(1,1)^3*V(1,2) + z(1,1)^3*M0*q;
    
else
% Equation to solve
b =zeros(1, N);
c = zeros(1, N-1);
b(1,1:N) =[ 1/h,2/h*ones(1,N/2-2),1/L/2+1/h,1/L/2+1/h,2/h*ones(1,N/2-2),1/h];
c(1,1:N-1) = [-1/h*ones(1,N/2-1),-1/2/L,-1/h*ones(1,N/2-1)];%c(N)=0 not used
a = c;
M =  sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c,a],N,N);
dydt = alpha*( -M*q);
end
end