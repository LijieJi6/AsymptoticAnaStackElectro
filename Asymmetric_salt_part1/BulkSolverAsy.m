function [temphi, temc2, temc1, j_curr, A, phi_uniform,phi_newL,phi_newR] = BulkSolverAsy( V, z, zeta,zetaR, N_all, len, index_fdm, x, h, x_ed,varphi_newL,varphi_newR,L)
% [temphi, temc2, temc1] = PoissonSolverAsy( z, eps, V, N_all, len, index_fdm, varphi_old, err_tol,x, x_ed, deltax)

temphi = zeros(1, N_all);
phi_uniform = zeros(1, N_all);
phi_newL  = zeros(1, N_all);
phi_newR  = zeros(1, N_all);
z2 = z(1,2);
z1 = z(1,1);
temc1 = -z2*ones(1,N_all);
temc2 = z1*ones(1,N_all);
n = size(zeta,2);
for ik=1:floor(length(len))
    x_1 = x(1,index_fdm(ik,1:len(1,ik)));
    % y = (x_1+ x_ed(1,ik))/eps;
    % deltay = deltax/eps;
    if ik==floor(length(len)/2)+1
        j_curr = (zeta(1,ik) - zetaR(1,ik) -V(1,1) + V(1,2))/2/L;
        A = V(1,1) - j_curr*x_ed(1,ik)-zeta(1,ik);
        phi_uniform(1,index_fdm(ik,1:len(1,ik))) = j_curr*x_1 - A;
    elseif ik<floor(length(len)/2)+1
        j_curr = (zeta(1,ik) - zeta(1,ik+1))/h;
        A = V(1,1) - j_curr*x_ed(1,ik)-zeta(1,ik);
        phi_uniform(1,index_fdm(ik,1:len(1,ik))) = j_curr*(x_1 -x_ed(1,ik)-x_ed(1,ik+1))-A;
    else
        j_curr = (zetaR(1,ik-1) - zetaR(1,ik))/h;
        A = V(1,2) - j_curr*x_ed(1,ik)-zetaR(1,ik-1);
        phi_uniform(1,index_fdm(ik,1:len(1,ik))) = j_curr*(x_1 -x_ed(1,ik)-x_ed(1,ik+1))-A;
        
    end
    temphi(1,index_fdm(ik,1:len(1,ik))) = j_curr*x_1 + A;
    phi_newL(1,index_fdm(ik,1:len(1,ik))) = varphi_newL(1,index_fdm(ik,1:len(1,ik))) +(j_curr*x_ed(1,ik)+A);
    phi_newR(1,index_fdm(ik,1:len(1,ik))) = varphi_newR(1,index_fdm(ik,1:len(1,ik))) +(j_curr*x_ed(1,ik+1)+A);
    
end
