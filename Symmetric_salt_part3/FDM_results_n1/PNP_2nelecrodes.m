clear
clc
tic
%% FDM solver
%%
H = 0;% assume
L = 1;
eps = (H+L)/100;
h = 25*eps;
n_ed = (floor(H/h) + 1)*2;
x_ed = zeros(1, n_ed);
n_0 =floor(n_ed/2);
deltax = 1/2000;1/2000;
xl = -L;
xr = L;
N = (xr + H - xl + H)/deltax;% x_l =-1, x_r = +1
T = 150;60;
V = [-0.1,0.1];
z = [1,-1];
% deltat=deltax^2;
deltat= deltax*0.2;
lambda = eps*deltat/(deltax^2)/2;
x = zeros(1,N);% half nodes on [-1,1]
for j=1:N
    x(j) = -L - H + (j-1/2)*deltax;
end
x_ed(1, 1: H/h +1) = ((1: 1:H/h +1) - 1)*h-L-H;
x_ed(1, H/h +2: 2*(H/h +1)) = ((H/h +2: 2*(H/h +1)) - H/h -2)*h +L;
% .. pre-compute all the index which is used to solve Poisson equaiton in
% each small interval
len = zeros(1,n_ed-1);
index_fdm = zeros(n_ed-1,floor(2/deltax));
for ik =1:floor(H/h)
    len(1,ik) = floor(h/deltax);
    index_fdm(ik,1:len(1,ik)) = floor(h/deltax)*(ik-1)+1:1:floor(h/deltax)*ik;
end
len(1, floor(H/h)+1) = floor(2*L/deltax);
index_fdm(floor(H/h)+1,1:len(1,floor(H/h)+1)) =  floor(H/deltax) +1:floor(H/deltax) + floor(2*L/deltax);
for ik = floor(H/h) +2:n_ed-1
    len(1,ik) = floor(h/deltax);
    index_fdm(ik,1:len(1,ik)) = floor(H/deltax) + floor(2*L/deltax)+ (ik-floor(H/h)-2)*floor(h/deltax)+1:floor(H/deltax)...
        +floor(2*L/deltax)+(ik-floor(H/h)-1)*floor(h/deltax);
end
t_store = [0.5,1,2,5,10,15,20,40,60,70,80,90,100];
pp = zeros(length(t_store)+1,N);
nn = zeros(length(t_store)+1,N);
pphi = zeros(length(t_store)+1,N);
rho = zeros(length(t_store)+1,N);
QQ = zeros(1,length(t_store)+1);
% phi = V(2)*x;
% temphi = V(2)*x;
% .. the initial values of potential
phi(1,1:floor(H/deltax)) = ones(1,floor(H/deltax))*V(1);
phi(1,floor(H/deltax + 2*L/deltax)+1:floor(2*H/deltax + 2*L/deltax)) = ones(1,floor(H/deltax))*V(2);
phi(1,floor(H/deltax)+1:floor(H/deltax + 2*L/deltax)) = 2*V(2)*x(1,floor(H/deltax)+1:floor(H/deltax + 2*L/deltax));
temphi = phi;
temc1 = -z(2)*ones(1,N);
temc2 = z(1)*ones(1,N);
WW1 = zeros(1,N);
WW2 = zeros(1,N);
a = zeros(1,N);
b = zeros(1,N);
c = zeros(1,N);%
d = zeros(N,1);
BW1 = zeros(1,N-1);
temBW1 = zeros(1,N-1);
BW2 = zeros(1,N-1);
temBW2 = zeros(1,N-1);
BBW1 = zeros(1,N-1);
BBW2 = zeros(1,N-1);
n = 0;
n_store =0;
n_store = n_store+1;
pphi(n_store,1:N) = temphi;
pp(n_store,1:N) = temc1;
nn(n_store,1:N) = temc2;
rho(n_store,1:N) = z(1)*pp(n_store,1:N)+z(2)*nn(n_store,1:N);
QQ(1,n_store) = sum(rho(n_store,1:N))*deltax;
t_q_store =[0.01:0.01:0.99,(1:0.5:T)];
Q_total = zeros(length(t_q_store),n_ed);
c_salt = zeros(1,length(t_q_store));
ik_store = 0;
while n < T/deltat
    WW1 = (3*z(1)*temphi - z(1)*phi)/2;
    WW2 = (3*z(2)*temphi- z(2)*phi)/2;
    temBW1 = z(1)*(temphi(1:N-1)+temphi(2:N))/2;
    BW1 = z(1)*(phi(1:N-1)+phi(2:N))/2;
    temBW2 = z(2)*(temphi(1:N-1)+temphi(2:N))/2;
    BW2 = z(2)*(phi(1:N-1)+phi(2:N))/2;
    BBW1 = (3*temBW1-BW1)/2;
    BBW2 = (3*temBW2-BW2)/2;
    % ... iteration for positive ion concentration
    b(1,1) = 1+lambda*exp(-BBW1(1)+WW1(1));
    b(1,N) = 1+lambda*exp(-BBW1(N-1)+WW1(N));
    b(1,2:N-1) = 1+lambda*exp(-BBW1(2:N-1)+WW1(2:N-1))+lambda*exp(-BBW1(1:N-2)+WW1(2:N-1));
    c(1,1) =-lambda*exp(-BBW1(1)+WW1(2));
    c(1,2:N-1) =-lambda*exp(-BBW1(2:N-1)+WW1(3:N));
    a(1,2:N-1) =-lambda*exp(WW1(1:N-2)-BBW1(1:N-2));
    a(1,N) =-lambda*exp(-BBW1(N-1)+WW1(N-1));
    d(1,1) =(1-lambda*exp(-BBW1(1)+WW1(1)))*temc1(1)+lambda*exp(-BBW1(1)+WW1(2))*temc1(2);
    d(N,1) =(1-lambda*exp(-BBW1(N-1)+WW1(N)))*temc1(N)+lambda*exp(-BBW1(N-1)+WW1(N-1))*temc1(N-1);
    d(2:N-1,1) =(1-lambda*exp(-BBW1(1:N-2)+WW1(2:N-1))-lambda*exp(-BBW1(2:N-1)+WW1(2:N-1))).*temc1(2:N-1)+...
        (lambda*exp(-BBW1(2:N-1)+WW1(3:N))).*temc1(3:N)+(lambda*exp(WW1(1:N-2)-BBW1(1:N-2))).*temc1(1:N-2);
    CoeffMatrix = sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],...
        [b(1,:),c(1,1:N-1),a(1,2:N)],N,N);
    temc1 = (CoeffMatrix\d)';% n+1 the time nodes
    % ... iteration for negative ion concentration
    b(1,1) = 1+lambda*exp(-BBW2(1)+WW2(1));
    b(1,N) = 1+lambda*exp(-BBW2(N-1)+WW2(N));
    b(1,2:N-1) = 1+lambda*exp(-BBW2(2:N-1)+WW2(2:N-1))+lambda*exp(-BBW2(1:N-2)+WW2(2:N-1));
    c(1,1) = -lambda*exp(-BBW2(1)+WW2(2));
    c(1,2:N-1) = -lambda*exp(-BBW2(2:N-1)+WW2(3:N));
    a(1,2:N-1) = -lambda*exp(WW2(1:N-2)-BBW2(1:N-2));
    a(1,N) = -lambda*exp(-BBW2(N-1)+WW2(N-1));
    d(1,1) = (1-lambda*exp(-BBW2(1)+WW2(1)))*temc2(1)+lambda*exp(-BBW2(1)+WW2(2))*temc2(2);
    d(N,1) = (1-lambda*exp(-BBW2(N-1)+WW2(N)))*temc2(N)+lambda*exp(-BBW2(N-1)+WW2(N-1))*temc2(N-1);
    d(2:N-1,1) = (1-lambda*exp(-BBW2(1:N-2)+WW2(2:N-1))-lambda*exp(-BBW2(2:N-1)+WW2(2:N-1))).*temc2(2:N-1)+...
        (lambda*exp(-BBW2(2:N-1)+WW2(3:N))).*temc2(3:N)+(lambda*exp(WW2(1:N-2)-BBW2(1:N-2))).*temc2(1:N-2);
    CoeffMatrix = sparse([1:1:N,1:1:N-1,2:1:N],[1:1:N,2:1:N,1:1:N-1],...
        [b(1,:),c(1,1:N-1),a(1,2:N)],N,N);
    temc2 = (CoeffMatrix\d)'; % n+1 the time nodes
    %...... for potential phi.... solve phi in each small interval
    phi = temphi;
    temphi = PoissonSolver(temc2, temc1, z, eps, V, N, len, index_fdm, deltax);
    n = n+1;
    
    
    % calculate the total diffuse charge around each electrode
    if sum((n*deltat) == t_q_store)==1
        ik_store =ik_store+1;
        tt_real(1,ik_store:ik_store) = (n*deltat);
        c_salt(1,ik_store) = z(1,1)*temc1(1,2000) -z(1,2)*temc2(1,2000);
        rho0 = z(1)*temc1+z(2)*temc2;
        Q_total(ik_store,floor(n_ed/2)) = sum(rho0(1,1:floor(0.5*h/deltax)),2)*deltax;
        if floor(n_ed/2)-1>1
            for ik=2:floor(n_ed/2)-1
                Q_total(ik_store,ik) = sum(rho0(1,floor((n_0-ik-0.5)*h/deltax)+1:1:floor((n_0-ik-0.5+1)*h/deltax)),2)*deltax;
            end
            Q_total(ik_store,1) = sum(rho0(1,floor(H/deltax) +1:floor(H/deltax) + floor(1*L/deltax)),2)*deltax+sum(rho0(1,floor((n_0-1-0.5)*h/deltax)+1:1:floor((n_0-1)*h/deltax)),2)*deltax;
        else
            Q_total(ik_store,1) = sum(rho0(1,floor(H/deltax) +1:floor(H/deltax) + floor(1*L/deltax)),2)*deltax;
        end
    end
    
    
    
    if sum((n*deltat) == [0.5,1,2,5,10,15,20,40,60,70,80,90,100])==1
        n_store = n_store+1;
        pphi(n_store,1:N) = temphi;
        pp(n_store,1:N) = temc1;
        nn(n_store,1:N) = temc2;
        rho(n_store,1:N) = z(1)*pp(n_store,1:N)+z(2)*nn(n_store,1:N);
        QQ(1,n_store) = sum(rho(n_store,1:N))*deltax;
    end
end

cl='FDMdata_symm_phi01';
outfile1=[cl,'.mat'];
save(outfile1,'pp','pphi','nn','rho','QQ','x','Q_total','tt_real','c_salt','-v7.3');
