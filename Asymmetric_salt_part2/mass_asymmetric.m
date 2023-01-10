function M = mass_asymmetric(t,q,z,len)
% Mass matrix function
len0 = floor(len/2);
C11=(z(1,1)*z(1,2)*exp(-z(1,1)*q(1:len0)) -z(1,1)*z(1,2)*exp(-z(1,2)*q(1:len0)));
C12 = (-sqrt(2)*sqrt(-z(1,2)*exp(-z(1,1)*q(1:len0)) +z(1,1)*exp(-z(1,2)*q(1:len0))+z(1,2)-z(1,1)));
C21=(z(1,1)*z(1,2)*exp(-z(1,1)*q(len0+1:2*len0)) -z(1,1)*z(1,2)*exp(-z(1,2)*q(len0+1:2*len0)));
C22 = (sqrt(2)*sqrt(-z(1,2)*exp(-z(1,1)*q(len0+1:2*len0)) +z(1,1)*exp(-z(1,2)*q(len0+1:2*len0))+z(1,2)-z(1,1)));
C11(2:len0) =2*C11(2:len0);
C21(1:len0-1) = 2*C21(1:len0-1);
M_d =[C11./C12;C21./C22]; 
M = sparse(1:2*(len0),1:2*(len0),M_d);
end