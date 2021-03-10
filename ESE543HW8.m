%% Question 1 
clear all 
clc
% State space model
A = [0 1; 1 0];
B = [0; 1];

% 1.a checking system stability
[vecA, valA] = eig(A);
% eigenvalues of A matrix are 1 and -1. The positive one makes the system unstable

% 1.b show unstable subspace is in controllable subspace
Mc = [B A*B];
symMc = sym(Mc);
rsMc = colspace(symMc);

% eigenvector for eigenvalue -1 for A matrix
vecA2 = vecA(:,2);
% linear combination for range space of Mc includes [1; 1]
% this is equivalent to the span of the positive eigenvector [0.707; 0.707]
% this is means that unstable subspace is in the controllable subspace

% 1.c state feedback to produce closed loop eigenvalues at -1, -1/2
% system is already in controllable canonical form
Kx = [-3/2 3/2];
Kx = -Kx;
% checking
Acl = A-B*Kx;
eigAcl = eig(Acl);

% 1.d Design a full order observer having poles at -4 and -5
% use the observer feedback to produce closed loop eigenvalues at -1/2, -1, -4, -5.
C = [1 0];
L = [9; 21];
of = [A-B*Kx B*Kx; zeros(2) A-L*C]; 
eigOF = eig(of);

% 1.e using a first-order observer with pole at -6. 
% Give a block diagram showing the controller as a single transfer function.
L1 = [1; 0];
L2 = [0; 1];
CP = [0 1];
CPAL2 = CP*A*L2;
CAL2 = C*A*L2;
K = [-6];
G1 = CPAL2 - K*CAL2;
G2 = CPAL2*K + CP*A*L1 - K*C*A*L1 - K*CAL2*K;
G3 = CP*B - K*C*B;
% Controller State Space
Ac = [G1-G3*Kx*L2];
Bc1 = [G2-G3*Kx*(L1+L2*K)];
Bc2 = [0];
Cc = [-Kx*L2];
Dc1 = [-Kx*(L1+L2*K)];
Dc2 = [0];
[num, denom] = ss2tf(Ac,Bc1,Cc,Dc1);

%% Question 2
clear all
clc
A = -1;
B = 1;
Q = 1;
R = 1;
N = 0;
[K, S, e] = lqr(A,B,Q,R,N);

syms P
eq = (P-sqrt(2)+1)/(P+sqrt(2)+1);
c = 0.1716;
RHS = c*exp(sqrt(2)*2);
Pt = simplify(solve(eq == RHS, P));
%% Question 3
clear all
clc

A = [0 1; 1 0];
B = [0; 1];
Q = [1 0; 0 0];
R = [1];
N = [0];
Mc = [B A*B];
rankMc = rank(Mc);
Qhalf = [1 0];
Mo = [Qhalf; Qhalf*A];
rankMo = rank(Mo);

syms l m n
P = [l m; m n];
riccati = P*A + transpose(A)*P + Q - P*B*inv(R)*transpose(B)*P;
riccati11 = solve(riccati(1,1) == 0, m);
[K, S, e] = lqr(A,B,Q,R,N);

l2 = 3.1076;
m2 = 2.4142;
n2 = 2.1974;

P = [l2 m2; m2 n2];
Kx = inv(R)*transpose(B)*P;

%% State Space
clear all 
clc
% plant state space
Ap = [0 1; 1 0];
Bp = [0; 1];
Cp = eye(2);
Dp = 0.*Cp*Bp;

% controller state space
R = [1];
l2 = 3.1076;
m2 = 2.4142;
n2 = 2.1974;
P = [l2 m2; m2 n2];
Kx = inv(R)*transpose(Bp)*P;
Ac = 0;
Bc1 = [0 0];
Bc2 = 0;
Cc = 0;
Dc1 = -Kx;
Dc2 = 0;

% Closed Loop State Space
%State Space for Closed Loop System
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl = Bcl(:,1);
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
Dcl = Dcl(:,1);
syscl = ss(Acl,Bcl,Ccl,Dcl);
% Lu
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

% Nyquist for Lu
figure(1);
nyquist(sys_u)
axis([-2 2 -2 2]);

% Bode for Lu
figure(2);
margin(sys_u)
freqResp = allmargin(sys_u);

% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));
figure(3)
sigma(sys_u,[],2);
title('Singular Values for I + Lu')
% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));
figure(4)
sigma(sys_u,[],3);
title('Singular Values for I + inv(Lu)')
% singular value stability margins for I+Lu
RDu_nGM = 1/(1+minsbar1);
RDu_pGM = 1/(1-minsbar1);
RDu_Pha = 2*asin(minsbar1/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi; 

% singular value stability margins for I+inv(Lu)
SDu_nGM = 1/(1+minsbar2);
SDu_pGM = 1/(1-minsbar2);
SDu_Pha = 2*asin(minsbar2/2);
SDu_nGM_dB = 20*log10(SDu_nGM);
SDu_pGM_dB = 20*log10(SDu_pGM);
SDu_Pha_deg = 180*RDu_Pha/pi; 
%% Riccati calculations
m1 = -0.4142;
m2 = 2.4142;
syms l n
eq31 = - n^2 + 2*m1;
eq32 = - n^2 + 2*m2;
n1 = solve(eq31 == 0, n);
n2 = solve(eq32 == 0, n);
n2 = 2.1974;
eq21 = l + n1 - m1*n1;
eq22 = l + n2 - m2*n2;
l2 = solve(eq22 == 0, l);
l2 = 3.1076;

