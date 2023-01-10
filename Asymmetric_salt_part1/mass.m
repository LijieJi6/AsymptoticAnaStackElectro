function M = mass(t,q,z,len)
% Mass matrix function
len0 = floor(length(len)/2)+1;
C1=sqrt(z(1,2)^2*z(1,1)-z(1,1)^2*z(1,2))./(1+(z(1,1)+z(1,2))/3*q);
C1(2:end-1) = C1(2:end-1)*2;
M_d =C1;
M = sparse(1:2*(len0),1:2*(len0),M_d);

end