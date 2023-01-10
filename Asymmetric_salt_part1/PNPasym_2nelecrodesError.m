clear
clc
tic
%% Asymptotic solver
%%
H = 0.5;% assume
L = 0.5;
eps = (H+L)/200;
h = 25*eps;
n_ed = (floor(H/h) + 1)*2;
x_ed = zeros(1, n_ed);
deltax = 1/2000;
xl = -L;
xr = L;
N = (xr + H - xl + H)/deltax;% x_l =-1, x_r = +1
T =150;
V = [-0.2,0.2];
t_store =0.1:0.1:T;% [0.5,1,2,5,10,15,20,40,60,70,80,90,100];
z = [2,-1];
n_0 =floor(n_ed/2);

% deltat=deltax^2;
deltat= deltax*0.2;
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
pp = zeros(length(t_store)+1,N);
nn = zeros(length(t_store)+1,N);
pphi = zeros(length(t_store)+1,N);
rho = zeros(length(t_store)+1,N);
QQ_L = zeros(1,length(t_store)+1);

% .. the initial values of potential
phi(1,1:floor(H/deltax)) = ones(1,floor(H/deltax))*V(1);
phi(1,floor(H/deltax + 2*L/deltax)+1:floor(2*H/deltax + 2*L/deltax)) = ones(1,floor(H/deltax))*V(2);
phi(1,floor(H/deltax)+1:floor(H/deltax + 2*L/deltax)) = 2*V(2)*x(1,floor(H/deltax)+1:floor(H/deltax + 2*L/deltax));
phi2(1,1:floor(H/deltax)) = ones(1,floor(H/deltax))*V(1);
phi2(1,floor(H/deltax + 2*L/deltax)+1:floor(2*H/deltax + 2*L/deltax)) = ones(1,floor(H/deltax))*V(2);
phi2(1,floor(H/deltax)+1:floor(H/deltax + 2*L/deltax)) = 2*V(2)*x(1,floor(H/deltax)+1:floor(H/deltax + 2*L/deltax));
temphi = phi;
temc1 = -z(1,2)*ones(1,N);
temc2 = z(1,1)*ones(1,N);


% only calculate the total diffuse charge in the left half domain
n = 0;
n_store =0;
n_store = n_store+1;
N_half = index_fdm(floor(length(len)/2)+1,len(1,floor(length(len)/2)+1))-floor(L/deltax);
pphi(n_store,1:N) = temphi(1,1:N);
pp(n_store,1:N) = temc1(1,1:N);
nn(n_store,1:N) = temc2(1,1:N);
rho(n_store,1:N) = z(1,1)*pp(n_store,1:N)+z(1,2)*nn(n_store,1:N);
% QQ_L(1,n_store) = sum(rho(n_store,1:N_half))*deltax;
t_q_store =[0.01:0.01:0.99,(1:0.5:T)];
Q_total = zeros(length(t_q_store),n_ed);
Q_totalR = zeros(length(t_q_store),n_ed);
ik_store = 0;

err_tol = power(10,-10);
zeta_asyR = zeros(1,floor(length(len)));
% .........ODE Solver ......................
%calculate [0,T_real] with linearized differential capacitance
T_real =0.05;T;
[tt1, zeta_ode1] = ODESolver(deltat,T,len, z, V, h,L,T_real);
% calculate [T_real,T]
% tt2=[];
% zeta_ode2=[];
y_initial =zeta_ode1(end,:);
[tt2, zeta_ode2] = ODESolver2(deltat,T,len, z, V, h,L,T_real,y_initial);
% zeta_ode =[zeta_ode1;zeta_ode2];
zeta_ode =[zeta_ode1;zeta_ode2];
tt = [tt1;tt2];
save('zeta_ode2000T150.mat','zeta_ode','tt','V','h')
% load zeta_ode
while n < T/deltat
    %........ solving all the (k,L) sub-domains, k=1, ..., n'
    zeta_asy= zeta_ode(n+2,1:end-1);% T/deltat rows, and 2 columns
    [varphi_newL, temc2_newL, temc1_newL] = PoissonSolverAsyL( z, eps, zeta_asy, N_half, len, index_fdm, x, err_tol, deltax);
    % ........ solving all the (k,R) sub-domains, k=1,2, ..., n
    zeta_asyR(1,1:floor(length(len))) = zeta_ode(n+2,2:end);
    [varphi_newR, temc2_newR, temc1_newR] = PoissonSolverAsyR( z, eps, zeta_asyR, N_half, len, index_fdm, x, err_tol, deltax);
    % ........solving the bulk areas ..........., k=1,2, ..., n
    [temphi_bulk, temc2_bulk, temc1_bulk, j_curr, A, phi_uniform,phi_newL,phi_newR] = BulkSolverAsy( V, z, zeta_asy,zeta_asyR, N_half, len, index_fdm, x, h, x_ed,varphi_newL,varphi_newR,L);
    
    n = n+1;
    %     n
    if sum((n*deltat) == t_store)==1
        n_store = n_store+1;
        pphi(n_store,1:N) = phi_newL + phi_newR + phi_uniform;
        pp(n_store,1:N) = temc1_newL + temc1_newR - abs(z(1,2))*ones(1,N);
        nn(n_store,1:N) = temc2_newL + temc2_newR- abs(z(1,1))*ones(1,N);
        rho(n_store,1:N) = z(1,1)*pp(n_store,1:N) + z(1,2)*nn(n_store,1:N);
        
    end
    
    % calculate the total diffuse charge around each electrode
    if sum((n*deltat) == t_q_store)==1
        ik_store =ik_store+1;
        tt_real(1,ik_store:ik_store) = (n*deltat);
        temc1_final = temc1_newL + temc1_newR - abs(z(1,2))*ones(1,N);
        temc2_final = temc2_newL + temc2_newR- abs(z(1,1))*ones(1,N);
        rho0 = z(1,1)*temc1_final +z(1,2)*temc2_final;
        Q_total(ik_store,floor(n_ed/2)) = sum(rho0(1,1:floor(0.5*h/deltax)),2)*deltax;
        for ik=2:floor(n_ed/2)-1
            Q_total(ik_store,ik) = sum(rho0(1,floor((n_0-ik-0.5)*h/deltax)+1:1:floor((n_0-ik-0.5+1)*h/deltax)),2)*deltax;
        end
        Q_total(ik_store,1) = sum(rho0(1,floor(H/deltax) +1:floor(H/deltax) + floor(1*L/deltax)),2)*deltax+sum(rho0(1,floor((n_0-1-0.5)*h/deltax)+1:1:floor((n_0-1)*h/deltax)),2)*deltax;
        
        Q_totalR(ik_store,floor(n_ed/2)) = sum(rho0(1,end-floor(0.5*h/deltax):end),2)*deltax;
        for ik=2:floor(n_ed/2)-1
            Q_totalR(ik_store,ik) = sum(rho0(1,floor((H+2*L)/deltax)+ floor((ik-1-0.5)*h/deltax)+1:1:floor((H+2*L)/deltax)+ floor((ik-0.5)*h/deltax)),2)*deltax;
        end
        Q_totalR(ik_store,1) = sum(rho0(1,floor((H+L)/deltax) +1:floor((H+L)/deltax) + floor((L+0.5*h)/deltax)),2)*deltax;
    end
    
end


% ...save data..............
cl='AsydataPNP_asymmError2000';
outfile1=[cl,'.mat'];
save(outfile1,'pp','pphi','nn','rho','x','Q_total','Q_totalR','tt_real','-v7.3');

toc