function M = mass1(t,q,z,len)
% Mass matrix function
len0 = floor(len/2);
C1=sqrt(z(1,2)^2*z(1,1)-z(1,1)^2*z(1,2))./(1+(z(1,1)+z(1,2))/3*q);
C1(2:end-1) = 2*C1(2:end-1); 
M = sparse(1:2*(len0),1:2*(len0),C1);
end