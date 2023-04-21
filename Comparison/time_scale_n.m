clear
clc
tic
% %% calculate the time scale at different H,L
mode0 = '3';
switch mode0
    case '1'
        H = 0.5;% assume
        L = 0.5;
    case '2'
        H = 2/3;
        L = 1/3;
    case '3'%'0.001'
        H = 1/1001;
        L =1000/1001;
end
ik = 0;
eig0 =zeros(1,125);
tau0 =zeros(1,125);
tau1 =zeros(1,125);
V =[-0.2,0.2];
z = [1,-1];
alpha = z(1,2)^2*z(1,1) - z(1,1)^2*z(1,2);
nn_ed0 = 2;
h0 = 0.5/(nn_ed0-1);
eps0 = h0/10;
T = 60;
deltax0 = eps0/10;
deltat= deltax0*0.2;
N = 2*(nn_ed0);
C_star0 = C_file(N,deltat,T, z, V, h0,L);
for nn_ed =2:2:250
    ik =ik+1;
    h = H/(nn_ed-1);
    n_ed = 2*(nn_ed);% total number of electrodes
    N = n_ed;
    % matrix to solve
    b =zeros(1, N);
    c = zeros(1, N-1);
    b(1,1:N) =[ 1/h,2/h*ones(1,N/2-2),1/L/2+1/h,1/L/2+1/h,2/h*ones(1,N/2-2),1/h];
    c(1,1:N-1) = [-1/h,-1/h*ones(1,N/2-2),-1/2/L,-1/h*ones(1,N/2-1)];%c(N)=0 not used
    a = [-1/h*ones(1,N/2-1),-1/2/L,-1/h*ones(1,N/2-2),-1/h];
    M =  sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c,a],N,N);
    M0 = M*h;
    C_star =ones(1,N)* C_star0(1,2);
    C_star(1,1) =  C_star0(1,1);
    C_star(1,end) =  C_star0(1,1);
    
    M_star = diag(C_star)\M0;
    eig_value2n = sort(eig(full(M_star)));
    tau0(1,ik) = h/(alpha*eig_value2n(2,1));
    tau1(1,ik) = h/(alpha*eig_value2n(1,1));
    
end
tau_com = sqrt(2)*(L+H)/L*tau0;
% ...save data..............
cl='TimeScaleNSymmetric';
outfile1=[cl,mode0,'.mat'];
save(outfile1,'H','L','tau0','tau1','tau_com','-v7.3');

toc


