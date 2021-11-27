%{
Source code for
{Y. Pan, X, Gao, X. Xu, Joint Estimation of Direction and Polarization for
Partiall Polarized Signals Using Tri-polarized Nested Array,
Radioengineering, Accept in November, 2021.}.
If you find this code useful, please consider citing the above paper.
%}
clear
source = deg2rad(linspace(-60,60,8)+90);
K = length(source);
lambda = 1;
dl = lambda/2;
L1 = 3;
L2 = 3;
pos = sort([(1:L1) (L1+1)*(1:L2)]).';
geoa = pos*dl;
M = length(geoa);
N = 500;
SNR = 10;
alpha = linspace(-pi/3,pi/3,K);  % POA
beta = linspace(-pi/5,pi/5,K);   % PEA
DOP = 0.9*ones(1,K); % DOP
sigma2c = DOP; % power of completely polarized, sigma2u+sigma2c = 1
sigma2u = 1 - sigma2c; % power of unpolarized

paras.K = K;
paras.M = M;
paras.N = N;
paras.dl = dl;
paras.lambda = lambda;
paras.L1 = L1;
paras.L2 = L2;
paras.pos = pos;
paras.geoa = geoa;

s = zeros(2*K,N);
B = zeros(M*3,2*K);
for k = 1:K
    G = [cos(alpha(k)) sin(alpha(k));-sin(alpha(k)) cos(alpha(k))];
    w = [cos(beta(k)); 1j*sin(beta(k))];
    Rs = sigma2u(k)/2*eye(2)+sigma2c(k)*G*(w*w')*G';
    [v, d] = eig(Rs);
    s((k-1)*2+1:2*k,:) = v*d.^(0.5)*(sqrt(1/2)*(randn(2,N)+1j*randn(2,N)));
    C =  [-1 0;0 sin(source(k));0 -cos(source(k))];
    B(:,(k-1)*2+1:2*k) = kron(exp(1j*2*pi*geoa/lambda*cos(source(k))),C);
end
noise = sqrt(1/(10^(SNR/10)))*sqrt(1/2)*(randn(3*M,N) + 1i*randn(3*M,N));
z = B*s+noise;

% eatimation of DOA, DOP, POA, and PEA
[ang,DOP_est,alpha_est,beta_est] = EfficientJDP(z,paras);
JDP_estimates = [ang,DOP_est,alpha_est,beta_est] 
