
function dydt = fODE(t, q, z, len, V, h,L)
N = 2*(floor(length(len)/2)+1);
len0 =(floor(length(len)/2)+1);
% Equation to solve
b =zeros(1, N);
c = zeros(1, N-1);
b(1,1:N) =[ 1/h,2/h*ones(1,len0-2),1/h+1/2/L,1/h+1/2/L,2/h*ones(1,len0-2),1/h];
c(1,1:N-1) = [-1/h*ones(1,len0-1),-1/2/L,-1/h*ones(1,len0-1)];%c(N)=0 not used
a = c;
M0 =  sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c,a],N,N);
dydt =(z(1,2)^2*z(1,1)-z(1,1)^2*z(1,2))*( [zeros(len0-1,1);-(V(1,2)-V(1,1))/2/L; -(V(1,1)-V(1,2))/2/L;...
    zeros(len0-1,1)] -M0*q);
end