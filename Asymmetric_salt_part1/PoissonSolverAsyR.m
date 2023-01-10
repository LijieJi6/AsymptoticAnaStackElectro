function [temphi, temc2, temc1] = PoissonSolverAsyR( z, eps, zeta, N_all, len, index_fdm, x, err_tol, deltax)
% [temphi, temc2, temc1] = PoissonSolverAsy( z, eps, V, N_all, len, index_fdm, varphi_old, err_tol,x, x_ed, deltax)
temphi = zeros(1, N_all);
varphi_old = zeros(1,N_all);
z2 = z(1,2);
z1 = z(1,1);
% n =size(zeta,2)-1;
for ik=1:floor(length(len))
    errPB = 10;
    deltay = deltax/eps;
    x_1 = x(1,index_fdm(ik,1:len(1,ik)));
    if ik == floor(length(len)/2)+1
        varphi_old(1,index_fdm(ik,1:len(1,ik))) = linspace(0, zeta(1,ik), length(x_1));
    else
        varphi_old(1,index_fdm(ik,1:len(1,ik))) = linspace(0, zeta(1,ik), length(x_1));
    end
    while errPB > err_tol
        N = len(1,ik);
        a = zeros(1,N);
        b =zeros(1, N);
        c = zeros(1, N);
        d = zeros(N,1);
        b(1,1) = 3 + deltay^2*(-z2*z1^2*exp(-z1*varphi_old(1,index_fdm(ik,1)))+z1*z2^2*exp(-z2*varphi_old(1,index_fdm(ik,1))));
        b(1,N) = 3 + deltay^2*(-z2*z1^2*exp(-z1*varphi_old(1,index_fdm(ik,len(1,ik))))+z1*z2^2*exp(-z2*varphi_old(1,index_fdm(ik,len(1,ik)))));
        b(1,2:N-1) = 2 + deltay^2*(-z2*z1^2*exp(-z1*varphi_old(1,index_fdm(ik,2:len(1,ik)-1)))+z1*z2^2*exp(-z2*varphi_old(1,index_fdm(ik,2:len(1,ik)-1))));
        c(1,1:N-1) = -1;%c(N)=0 not used
        a(1,2:N) = -1;%a(1ï¼?0 not used
        d(1:N,1) = deltay^2*(-z2*z1*(exp(-z1*varphi_old(1,index_fdm(ik,1:len(1,ik))))...
            +z1*varphi_old(1,index_fdm(ik,1:len(1,ik))).*exp(-z1*varphi_old(1,index_fdm(ik,1:len(1,ik)))))...
            +z1*z2*(exp(-z2*varphi_old(1,index_fdm(ik,1:len(1,ik))))+...
            z2*varphi_old(1,index_fdm(ik,1:len(1,ik))).*exp(-z2*varphi_old(1,index_fdm(ik,1:len(1,ik))))));
        if ik==floor(length(len)/2)+1
            d(N,1) = 2*zeta(1,ik)+ d(N,1);
            d(1,1) = 2*0+ d(1,1);
        else
            d(N,1) = 2*zeta(1,ik)+ d(N,1);
            d(1,1) = 2*0+ d(1,1);
        end
        CoeffMatrix = sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c(1,1:N-1),a(1,2:N)],N,N);
        temphi(1,index_fdm(ik,1:len(1,ik))) = (CoeffMatrix\d)';
        errPB = max(abs(temphi(1,index_fdm(ik,1:len(1,ik))) - varphi_old(1,index_fdm(ik,1:len(1,ik)))));
        varphi_old(1,index_fdm(ik,1:len(1,ik))) = temphi(1,index_fdm(ik,1:len(1,ik)));
    end
    
end
temc2 = z1*exp(-z2*temphi);
temc1 = -z2*exp(-z1*temphi);