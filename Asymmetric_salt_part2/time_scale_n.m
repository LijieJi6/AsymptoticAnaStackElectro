clear
clc
tic
% %% calculate the time scale at different H,L
H = 0.5;% assume
L = 0.5;
ik = 0;
tau0 =zeros(1,125);
tau1 =zeros(1,125);
V =[-0.2,0.2];
z = [2,-1];
alpha = z(1,2)^2*z(1,1) - z(1,1)^2*z(1,2);
nn_ed0 =5;
h = H/(nn_ed0-1);
eps = h/10;
deltax = eps/10;
deltat= deltax*0.2;
T_real = 0.5;
T =10;
n_ed0 = 2*(nn_ed0);% total number of electrodes
N = n_ed0;
[q0,C_star0] = C_file_asymmetric(N,deltat,T, z, V, h,L,T_real);

mode0 = '1';
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


for nn_ed = 2:2:250
    ik =ik+1;
    h = H/(nn_ed-1);
    n_ed = 2*(nn_ed);% total number of electrodes
    % eps = h/25;
    % deltax = eps/10;
    % deltat= deltax*0.2;
    N = n_ed;
    C_star =ones(1,N);
    % matrix to solve
    b =zeros(1, N);
    c = zeros(1, N-1);
    b(1,1:N) =[ 1/h,2/h*ones(1,N/2-2),1/L/2+1/h,1/L/2+1/h,2/h*ones(1,N/2-2),1/h];
    c(1,1:N-1) = [-1/h,-1/h*ones(1,N/2-2),-1/2/L,-1/h*ones(1,N/2-1)];%c(N)=0 not used
    a = [-1/h*ones(1,N/2-1),-1/2/L,-1/h*ones(1,N/2-2),-1/h];
    M =  sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],[b(1,:),c,a],N,N);
    M0 = M*h;
    % the differential capacitance should be calculated in two parts
    
    % T_real = 0.5;
    % [q0,C_star0] = C_file_asymmetric(N,deltat,T, z, V, h,L,T_real);
    
    C_star(1,1:nn_ed) =ones(1,nn_ed)* C_star0(1,2);
    C_star(1,1) =  C_star0(1,1);
    C_star(1,nn_ed+1:end) =  ones(1,nn_ed)*C_star0(1,end-1);
    C_star(1,end) = C_star0(1,end);
    
    M_star = diag(C_star)\M0;
    eig_value2n = sort(eig(full(M_star)));
    tau0(1,ik) = h/(alpha*eig_value2n(2,1));
    tau1(1,ik) = h/(alpha*eig_value2n(1,1));
end
% ...save data..............
cl='TimeScaleNASymmetric';
outfile1=[cl,mode0,'.mat'];
save(outfile1,'H','L','tau0','-v7.3');
plot(2:2:250,tau0)
toc
