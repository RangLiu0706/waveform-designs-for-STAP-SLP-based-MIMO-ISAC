% The proposed MM-neADMM algorithm for constant-modulus waveform design.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint waveform and filter designs for STAP-SLP-based MIMO-DFRC systems,”IEEE J. Sel. Areas Commun., vol. 40, no. 6, pp. 1918-1931, Jun. 2022.
% Download this paper at: https://ieeexplore.ieee.org/document/9769997
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-01-31
% Inputs: Prms: the structure of system parameters;
%         H0: the channel;
%         S: the transmitted communication symbols;
%         rho: the penalty;
% Outputs: x: transmit waveform; 
%          SINRiter: the achieved radar SINR
function [x,SINRiter] = get_x_my_CM(Prms,H0,S,rho)

Ku = Prms.Ku; Nt = Prms.Nt; Nr = Prms.Nr; M = Prms.M; N = Prms.N; Phi = Prms.Phi;
Nmax = Prms.Nmax; res_th = Prms.res_th; sigma2 = Prms.sigma2;
A = Prms.A; P = Prms.P; beta = Prms.beta;

%%% construct H
H = zeros(2*Ku*M*N,M*N*Nt);
for n = 1:1:M*N
    en = zeros(1,M*N);
    en(n) = 1;
    for ku = 1:1:Ku
        t1 = exp(-1j*angle(S(ku,n)))*( sin(Phi)+exp(-1j*pi/2)*cos(Phi) );
        t2 = exp(-1j*angle(S(ku,n)))*( sin(Phi)-exp(-1j*pi/2)*cos(Phi) );
        H((2*ku-2)*M*N+n,:) = kron(en,t1.*H0(ku,:));
        H((2*ku-1)*M*N+n,:) = kron(en,t2.*H0(ku,:));
    end
end

%%% initialize x
cvx_begin quiet
variable x(M*N*Nt,1) complex
variable t
maximize t
subject to
t <= real(H*x);
abs(x) <= sqrt(P/(M*N*Nt));
cvx_end

X = x*x';
Phi_Siter = sum( pagemtimes( pagemtimes(A(:,:,1:end-1),X),'none',A(:,:,1:end-1),'ctranspose'), 3);
T1 = A(:,:,end)'/(Phi_Siter+sigma2*eye(M*N*Nr));
bt = T1*A(:,:,end)*x;
T2 = T1'*X*T1;
Dt = sum( pagemtimes( pagemtimes(A(:,:,1:end-1),'ctranspose',T2,'none'), A(:,:,1:end-1)), 3);

iter = 0;
res = 1;
vobj = zeros(1,Nmax);
vres = zeros(1,Nmax);
SINRiter = zeros(1,Nmax);
mu = zeros(M*N*Nt,1);
lambda = zeros(M*N*Nt,1);
y = x;

while iter < Nmax && res > res_th

    TT = Dt + rho/2*eye(M*N*Nt);
    for i = 1:1:M*N*Nt
        TT(i,i) = real(TT(i,i));
    end
    RR = chol(TT);

    cvx_begin quiet
    variable x(M*N*Nt,1) complex
    minimize real(x'*(RR'*RR)*x)-real((bt'+rho*y'-lambda')*x)
    subject to
    real(H*x) >= beta;
    abs(x) <= sqrt(P/(M*N*Nt));
    cvx_end

    X = x*x';
    Phi_Siter = sum( pagemtimes( pagemtimes(A(:,:,1:end-1),X),'none',A(:,:,1:end-1),'ctranspose'), 3);
    T1 = A(:,:,end)'/(Phi_Siter+sigma2*eye(M*N*Nr));
    bt = T1*A(:,:,end)*x;
    T2 = T1'*X*T1;
    Dt = sum( pagemtimes( pagemtimes(A(:,:,1:end-1),'ctranspose',T2,'none'), A(:,:,1:end-1)), 3);

    %%% update y
    y =  0.5*( abs(x+lambda/rho) + real( sqrt(P/(M*N*Nt))-mu/rho ) ).*exp(1j*angle(x+lambda/rho));
    %%% update mu, lambda
    mu = mu + rho*(abs(y)-sqrt(P/(M*N*Nt)));
    lambda = lambda + rho*(x-y);

    iter = iter + 1;
    vres(iter) = norm(x-y,2)^2 + norm(abs(y)-sqrt(P/(M*N*Nt)),2)^2;
    vobj(iter) = real(-x'*bt) + ( norm(abs(y)-sqrt(P/(M*N*Nt))+mu/rho,2)^2 + norm(x-y+lambda/rho,2)^2 )*rho/2;
    SINRiter(iter) = 10*log10(real(sigma2*x'*bt));
    if iter > 1
        res = max([abs(1-SINRiter(iter)/SINRiter(iter-1)),abs(1-vobj(iter)/vobj(iter-1))] );
    end
end
vobj(iter+1:end) = [];
vres(iter+1:end) = [];
SINRiter(iter+1:end) = [];

