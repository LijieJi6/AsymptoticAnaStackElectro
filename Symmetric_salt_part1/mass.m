function M = mass(t,q,z,len)
% Mass matrix function
M_d = sqrt(2)*power(abs(z(1,1)),3/2)*cosh(z(1,1)*q/2);
M_d(1:end-1)=2*M_d(1:end-1);
M = sparse(1:floor(length(len)/2)+1,1:floor(length(len)/2)+1,M_d);

end