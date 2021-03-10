%%
clear all
clc
num = -1;
denom = [1 1];
[A,B,C,D] = tf2ss(num, denom);
sys1 = tf2ss(num, denom);
tf1 = tf(num,denom);
figure(1)
margin(tf1)
sv = sigma(tf1,[],2);
smallgain = 1./sv;
figure(2)
plot(sv)
figure(3)
plot(smallgain)


%%
clear all
clc
num = 1000;
dem = 1;
tf1 = tf(num, dem);
sigma(tf1,[],2);

%%
clear all
clc

syms s A B F L C N M
I = eye(2);
Amat = [A-(B*F) B*F; 0 A-(L*C)];
Bmat = [0; (L*N)+M];
Cmat = [0 1];
TF1 = Cmat*inv((s.*I) - Amat)*Bmat; 
ilaplace(TF1)

%% Question 1 d-e
clear all

% Initiating Constants
ZaV = -1.21;
ZdV = -0.1987;
Ma = 44.2506;
Md = -97.2313;
% Creating plant state space
Ap = [ZaV 1; Ma 0];
Bp = [ZdV; Md];
% Tracking the angle of attack
Cp = [1 0];
Dp = 0.*Cp*Bp;

% Creating wiggle system
w0 = 1;
Aw = [ 0 1 0 0
-w0*w0 0 Cp
0.*ones(2,2) Ap];
Bw = [ 0 0 Bp']';

% Checking if Aw and Bw is controllable
Kalman = [Bw Aw*Bw Aw^2*Bw Aw^3*Bw];
rankKal = rank(Kalman);
% Rank of Kalman matrix is 4 which is same size as row of Aw so system is
% controllable

% Placing Negative Poles
% Desired poles
dp1 = -5;
dp2 = -6;
dp3 = -7;
dp4 = -8;
% Placing poles
Kp = place(Aw, Bw, [dp1, dp2, dp3, dp4]);
% Creating new wiggle matrix with desired poles
nAw = Aw -(Bw*Kp);
eigNAw = eig(nAw);
% Eigenvalues are the desired poles so pole placement was done correctly

% Creating linear quadratic regulator(LQR)
R = 1;
Q = 0.*nAw;
Q(1,1) = 100000.;
Q(2,2) = 10000.;
[K,S,E]=lqr(nAw,Bw,Q,0.001*R);

Acl = nAw-Bw*K;
eig(Acl)