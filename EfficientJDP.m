%{
Source code for
{Y. Pan, X, Gao, X. XU, Joint Estimation of Direction and Polarization for
Partiall Polarized Signals Using Tri-polarized Nested Array,
Radioengineering, Accept in November, 2021.}.
If you find this code useful, please consider citing the above paper.
%}
function [ang,DOP_est,alpha_est,beta_est] = EfficientJDP(z,paras)
K = paras.K ;
M = paras.M ;
L1 = paras.L1 ;
L2 = paras.L2 ;
lambda = paras.lambda ;
dl = paras.dl;
N = paras.N;
pos = paras.pos ;
geoa = paras.geoa;

x1 = z(1:3:end,:); x2 = z(2:3:end,:); x3 = z(3:3:end,:);
Rxx = x1*x1'/N; Rxy = x1*x2'/N; Ryz = x2*x3'/N; Ryy = x2*x2'/N; Rxz = x1*x3'/N; Rzz = x3*x3'/N;
R = Rxx + Ryy + Rzz; % Sub-covariance Addition
ind = pos - pos.';
uC = uniquetol(ind);
Mp = length(uC);
vx = zeros(Mp,1);
for i = 1:Mp
    vx(i) = mean(R(ind==uC(i))); % outputs of the virtual ULA
end
midVx = (Mp+1)/2;
Tvxidx = (midVx-1:-1:0) + (1:midVx).';
Tvx = vx(Tvxidx); % Toeplitz Matrix for Rank Restoration 

% DOA estimation
[V,D] = eig(Tvx);
[~, order] = sort(abs(diag(D)));
Un = V(:,order(1:midVx-K));
C = Un*Un';
ind = (midVx:-1:1).' + (0:1:midVx-1);
a = accumarray(ind(:),conj(C(:)));
ra = roots(a);
rb = ra(abs(ra)<1);
[~,I] = sort(abs(abs(rb)-1));
ang = sort(acosd(angle(rb(I(1:K)))*lambda/(2*pi*dl)));

% Polarization estimation
Rz = z*z'/N;
[~,De] = eig(Rz);
[n,~] = sort(abs(diag(De)));
esNoise = mean(n(1:2*(L1+L2)-K));  % estimation of noise power
A1 = zeros(M^2,K);
A2 = zeros(3*(M^2),K);
for k = 1:K
    a = exp(1j*2*pi*geoa/lambda*cosd(ang(k)));
    A1tmp = kron(conj(a),a);
    A1(:,k) = A1tmp;
    A2(:,k) = kron([sind(ang(k))^2,cosd(ang(k))^2,-sind(ang(k))*cosd(ang(k))].',A1tmp);
end
p11 = A1\vec(Rxx-esNoise*eye(M));
p22 = A2\[vec(Ryy - esNoise*eye(M));vec(Rzz - esNoise*eye(M));vec(Ryz)];
Phi = diag(exp(1j*ang/180*pi));
p12 = (A1*Phi)\(vec(Rxz)-vec(1j*Rxy));
P = [p11 conj(p12) p12 p22];
pu = zeros(K,1); pc = zeros(K,1); omega = zeros(K,1); 
for k = 1:K
    [v, d] = eig(reshape(P(k,:),2,2));
    [dd, ind] = sort(abs(diag(d)));
    pu(k) = dd(1)*2;
    pc(k) = dd(2)-dd(1);
    e = v(:,ind(1));
    omega(k) = (1j*e(2)-e(1))/(1j*e(2)+e(1));
end
DOP_est = pc./(pu+pc);
alpha_est = angle(omega)/2;
beta_est = pi/4-atan(abs(omega));