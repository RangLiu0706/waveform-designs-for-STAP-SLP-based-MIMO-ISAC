% This Matlab script can be used to generate Fig. 5 in the paper:
% R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint waveform and filter designs for STAP-SLP-based MIMO-DFRC systems,”IEEE J. Sel. Areas Commun., vol. 40, no. 6, pp. 1918-1931, Jun. 2022.
% Download this paper at: https://ieeexplore.ieee.org/document/9769997
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-31
clc
clear

Nt = 6; %%% the number of transmit antennas
Nr = 6; %%% the number of receive antennas
L = 2; %%% five clutter range cells
sigma2 = 1; %%% noise power
sigman2 = 0.01;  %%% communication noise power
sigma_c2 = 1; %%% clutter power
sigma_t2 = 1; %%% target power
M = 4; %%% the number of pulses
N = 8; %%% the length of the waveform
Jl_hat = zeros(N,N,L); %%% shift matrix  -L---L
Jl = zeros(M*Nr*N,M*Nr*N,L);
for l = 1:1:L
    for i = 1:1:N
        for j = 1:1:N
            if i-j+l == 0
                Jl_hat(i,j,l) = 1;
            end
        end
    end
    Jl(:,:,l) = kron(eye(M*Nr),Jl_hat(:,:,l)');
end
T = zeros(N*Nt,N*Nt);  %%% permutation matrix
for i = 1:1:N*Nt
    m = floor(i/N);
    n = i - m*N;
    if n == 0
        temp =  (N-1)*Nt + m;
    else
        temp = (n-1)*Nt + m + 1;
    end
    T(i,temp) = 1;
end
%%% clutter
Nc = 60; %%% the numebr of clutter patches in each range cell
Mcl = zeros(M*Nr*Nt,M*Nr*Nt); %%% inner covariance matrix
for ctemp = 1:1:Nc
    theta_cl = ctemp/Nc*2*pi;
    rcl = exp(1i*pi*sin(theta_cl)*(0:1:Nr-1).');
    tcl = exp(1i*4*pi*sin(theta_cl)*(0:1:Nt-1).');
    dcl = exp(1i*pi*sin(theta_cl)*(0:1:M-1).');
    ucl = kron(kron(dcl,rcl),tcl);
    Mcl = Mcl + sigma_c2*(ucl*ucl');
end
rlr = rank(Mcl);
[Vcl,Dcl] = eigsort(Mcl);
A = zeros(M*Nr*N,M*Nt*N,rlr*(2*L+1)+1);
for r = 1:1:rlr
    ur = sqrt(Dcl(r,r))*Vcl(:,r);
    Ar = zeros(M*N*Nr,M*N*Nt);
    for m = 1:1:M
        Urm = reshape( ur((m-1)*Nr*Nt+1:m*Nr*Nt), Nt, Nr);
        Ar((m-1)*N*Nr+1:m*N*Nr,(m-1)*N*Nt+1:m*N*Nt) = kron(Urm.',eye(N))*T;
    end
    for l = -L:1:-1
        A(:,:,(r-1)*(2*L+1)+L+1+l) = Jl(:,:,-l)'*Ar;
    end
    A(:,:,(r-1)*(2*L+1)+L+1) = Ar;
    for l = 1:1:L
        A(:,:,(r-1)*(2*L+1)+L+1+l) = Jl(:,:,l)*Ar;
    end
end
%%% target
theta = 0;  %%% target DoA
fd = 0.3;   %%% target Doppler frequency
r0 = exp(1i*pi*sin(theta)*(0:1:Nr-1).'); %%% receive steering vector
t0 = exp(1i*4*pi*sin(theta)*(0:1:Nt-1).'); %%% transmit steering vector
d0 = exp(1i*2*pi*fd*(0:1:M-1).');  %%% the target Doppler response vector
u0 = kron(kron(d0,r0),t0); %%% spatial-temporal steering vector
for m = 1:1:M
    U0m = reshape( u0((m-1)*Nr*Nt+1:m*Nr*Nt), Nt, Nr);
    A((m-1)*N*Nr+1:m*N*Nr,(m-1)*N*Nt+1:m*N*Nt,end) = kron(U0m.',eye(N))*T;
end
Phi = pi/2;  %%% BPSK modulated
Pt = 100;  %%% total transmit power
Nmax = 500;   %%%% maximum iterations
res_th = 1e-4;%%%% convergence tolerance
SNR = 5;

Ku_range = (1:1:6);
N_sim = 100;
SINR_my_CE = zeros(1,length(Ku_range));
SINR_my_PAR = zeros(1,length(Ku_range));
SINR_my_CES = zeros(1,length(Ku_range));

epsi = 1.8*sqrt(Pt/M/N/Nt); %%% similarity level
ratio = 2; %%% PAPR level
Prms.Nt = Nt; Prms.Nr = Nr; Prms.M = M; Prms.N = N;
Prms.Phi = Phi; Prms.Nmax = Nmax; Prms.res_th = res_th; Prms.sigma2 = sigma2;
Prms.sigma_t2 = sigma_t2; Prms.A = A; Prms.P = Pt; Prms.epsi = epsi;
Prms.beta = sin(Phi)*sqrt(sigman2*10^(0.1*SNR)); Prms.ratio = ratio;
rho = 1;%%% penalty
load('H_K.mat')
for sim = 1:1:N_sim
    tic
    sim
    %     H0 = sqrt(0.5)* (rand(Ku_range(end),Nt) + 1j*rand(Ku_range(end),Nt) );
    %     S0 = exp(1i*(Phi+2*Phi*randi([0,pi/Phi-1],Ku_range(end),M*N)));
    H0 = HH(:,:,sim);
    S0 = SS(:,:,sim);

    for Ku_index = 1:1:length(Ku_range)
        Ku = Ku_range(Ku_index)
        Prms.Ku = Ku;

        H = H0(1:Ku,:);
        S = S0(1:Ku,:);

        [x_my_ce,SINRiter_my_ce] = get_x_my_CM(Prms,H,S,rho);
        SINR_my_CE(Ku_index) = SINR_my_CE(Ku_index) + max(SINRiter_my_ce);

        [x_my_ces,SINRiter_my_ces] = get_x_my_CMS(Prms,H,S,rho*5);
        SINR_my_CES(Ku_index) = SINR_my_CES(Ku_index) + max(SINRiter_my_ces);

        [x_my_par,SINRiter_my_par] = get_x_my_PAPR(Prms,H,S,1);
        SINR_my_PAR(Ku_index) = SINR_my_PAR(Ku_index) + max(SINRiter_my_par);

    end
    toc
end

SINR_my_CE = SINR_my_CE/sim;
SINR_my_PAR = SINR_my_PAR/sim;
SINR_my_CES = SINR_my_CES/sim;

figure;
plot(Ku_range,SINR_my_CE,'-o','color',[0.7,0,0],'LineWidth',1.5)
hold on
plot(Ku_range,SINR_my_CES,'-s','color',[0,0.5,0],'LineWidth',1.5)
plot(Ku_range,SINR_my_PAR,'-^','color',[0,0,0.7],'LineWidth',1.5)
hold off
xlabel('Number of users {\it K}_u');
ylabel('Radar SINR (dB)');
grid on
legend('Proposed, CM','Proposed, CMS','Proposed, PAPR')
