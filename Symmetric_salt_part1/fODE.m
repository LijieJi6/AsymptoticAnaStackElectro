
function dydt = fODE(t, q, z, len, V, h,L)
N = floor(length(len)/2)+1;
if N==1
    M0 = 1;
    dydt = z(1,1)^3*V(1,2) + z(1,1)^3*M0*q;
    
else
    % Equation to solve
    b =zeros(1, N);
    c = zeros(1, N-1);
    b(1,1:N) =[ 1/L+1/h,2/h*ones(1,N-2),1/h];
    c(1,1:N-1) = -1/h;%c(N)=0 not used
    a = c;
    M0 =  sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c,a],N,N);
    dydt = [-2*z(1,1)^3*V(1,2)/L;zeros(floor(length(len)/2),1)] - 2*z(1,1)^3*M0*q;
end
end