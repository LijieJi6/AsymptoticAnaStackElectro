function temphi = PoissonSolver(temc2, temc1, z, eps, V, N_all, len, index_fdm, deltax)

temphi = zeros(1, N_all);
for ik=1:length(len)
    N = len(1,ik);
    a = zeros(1,N);
    b =zeros(1, N);
    c = zeros(1, N);
    d = zeros(N,1);
    b(1,1) = 3*eps^2;
    b(1,N) = 3*eps^2;
    b(1,2:N-1) = 2*eps^2;
    c(1,1:N-1) = -eps^2;%c(N)=0 not used
    a(1,2:N) = -eps^2;%a(1）=0 not used
    d(2:N-1,1) = deltax^2*(z(1)*temc1(1,index_fdm(ik,2:len(1,ik)-1))+z(2)*temc2(1,index_fdm(ik,2:len(1,ik)-1)));
    if ik <= floor(length(len)/2)
        d(1,1) = 2*eps^2*V(1)+deltax^2*(z(1)*temc1(1,index_fdm(ik,1))+z(2)*temc2(1,index_fdm(ik,1)));
        d(N,1) = 2*eps^2*V(1)+deltax^2*(z(1)*temc1(1,len(1,ik))+z(2)*temc2(1,len(1,ik)));
    elseif ik == floor(length(len)/2)+1
        d(1,1) = 2*eps^2*V(1)+deltax^2*(z(1)*temc1(1,index_fdm(ik,1))+z(2)*temc2(1,index_fdm(ik,1)));
        d(N,1) = 2*eps^2*V(2)+deltax^2*(z(1)*temc1(1,len(1,ik))+z(2)*temc2(1,len(1,ik)));
    else
        d(1,1) = 2*eps^2*V(2)+deltax^2*(z(1)*temc1(1,index_fdm(ik,1))+z(2)*temc2(1,index_fdm(ik,1)));
        d(N,1) = 2*eps^2*V(2)+deltax^2*(z(1)*temc1(1,len(1,ik))+z(2)*temc2(1,len(1,ik)));
    end
    
    CoeffMatrix = sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c(1,1:N-1),a(1,2:N)],N,N);
    temphi(1,index_fdm(ik,1:len(1,ik))) = (CoeffMatrix\d)';
    
end