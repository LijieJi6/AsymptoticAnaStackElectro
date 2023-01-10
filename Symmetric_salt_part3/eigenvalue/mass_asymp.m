function M = mass_asymp(t,q,z,len)
% Mass matrix function
len0 = floor(len/2);

M_d = sqrt(2)*power(abs(z(1,1)),3/2)*cosh(z(1,1)*q/2);
M_d(2:end-1) = 2*M_d(2:end-1);
M = sparse(1:2*(len0),1:2*(len0),M_d);

end