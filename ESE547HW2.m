%% Question 1
clear all
clc
%*************************************************************************
% Plant Model
%*************************************************************************
Ap = [ -0.576007 -3255.07 4.88557 9.25796;
-0.0410072 -0.488642 -2.03681 0 ;
0 0 0 1 ;
0 0 -8882.64 -133.266 ];
Bp = [ 0 ; 0 ; 0 ; 8882.64];
Cp = [ 1 0 0 0;
0 1 0 0];
Dp = 0.*Cp*Bp;
%*******************************************************
% Static Output Feedback With New Form To Test Prior To Substituting The Big Plant Model
%*******************************************************
%Close the loop to test the model
% Plant form xdot = Apx + Bpu; y = Cpx +Dpu
% Controller xcdot = Acxc + Bc1y + Bc2r
% u = Ccxc + Dc1y + Dc2r
Ac = [ 0 ];
Bc1 = [ -1 0];
Bc2 = [ 1 ];
Cc = [ 0.0107349];
Dc1 = [ -0.0411729 11.4003];
Dc2 = [ 0 ];

%Closed loop state space
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl = Bcl(:,1);
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
Dcl = Dcl(:,1);
syscl = ss(Acl,Bcl,Ccl,Dcl);

%SS model of loop gain Lu at the plant input
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

%SS model of loop gain Ly at the plant output
Aout= [Ap Bp*Cc; 0.*Bc1*Cp Ac];
Bout = [Bp*Dc1; Bc1];
Cout= -[Cp Dp*Cc];%change sign for loop gain
Dout= -[Dp*Dc1];
sys_y= ss(Aout,Bout,Cout,Dout);

% Nyquist for Lu
figure(1);
nyquist(sys_u)
axis([-2.5 2 -2 2]);

% Bode for Lu
figure(2);
margin(sys_u)
freqResp = allmargin(sys_u);

% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));
figure(3);
sigma(sys_u,[],2);
% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));
figure(4);
sigma(sys_u,[],3);
% singular value stability margins
RDu_nGM1 = 1/(1+minsbar1);
RDu_pGM1 = 1/(1-minsbar1);
RDu_Pha1 = 2*asin(minsbar1/2);
RDu_nGM_dB1 = 20*log10(RDu_nGM1);
RDu_pGM_dB1 = 20*log10(RDu_pGM1);
RDu_Pha_deg1 = 180*RDu_Pha1/pi; 

RDu_nGM2 = 1/(1+minsbar2);
RDu_pGM2 = 1/(1-minsbar2);
RDu_Pha2 = 2*asin(minsbar2/2);
RDu_nGM_dB2 = 20*log10(RDu_nGM2);
RDu_pGM_dB2 = 20*log10(RDu_pGM2);
RDu_Pha_deg2 = 180*RDu_Pha2/pi; 

%% Question 3
clear all
clc

% Constants
Ka = -0.0015;
Kq = -0.32;
az = 2.0;
aq = 6.0;
V = 886.78; %fps
ZaV = -1.3046;
ZsV = -0.2142;
Ma = 47.7109;
Ms = -104.8346;
Za = ZaV*V;
Zs = ZsV*V;

% plant state space
Ap = [ZaV 1; Ma 0];
Bp = [ZsV; Ms];
Cp = [Za 0; 0 1];
Dp = [Zs; 0]; 

%State Space for Controller
Ac = [0 0; Kq*aq 0];
Bc1 = [-Ka*az 0; -Ka*Kq*aq -Kq*aq];
Bc2 = [Ka*az; Ka*Kq*aq];
Cc = [Kq 1];
Dc1 = [-Ka*Kq -Kq];
Dc2 = [Ka*Kq];

%State Space for Closed Loop System
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl = Bcl(:,1);
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
Dcl = Dcl(:,1);
syscl = ss(Acl,Bcl,Ccl,Dcl);

%Step response plot
figure(1);
step(syscl)

%SS model of loop gain Lu at the plant input
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

% Lu stability margin
% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));

% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));

% singular value stability margins
RDu_nGM1 = 1/(1+minsbar1);
RDu_pGM1 = 1/(1-minsbar1);
RDu_Pha1 = 2*asin(minsbar1/2);
RDu_nGM_dB1 = 20*log10(RDu_nGM1);
RDu_pGM_dB1 = 20*log10(RDu_pGM1);
RDu_Pha_deg1 = 180*RDu_Pha1/pi; 

RDu_nGM2 = 1/(1+minsbar2);
RDu_pGM2 = 1/(1-minsbar2);
RDu_Pha2 = 2*asin(minsbar2/2);
RDu_nGM_dB2 = 20*log10(RDu_nGM2);
RDu_pGM_dB2 = 20*log10(RDu_pGM2);
RDu_Pha_deg2 = 180*RDu_Pha2/pi; 

%% Actuator dynamics
% syms s t
% s = tf('s');
% actuatorTF = 1/((t*s) + 1);
% delta = simplify(actuatorTF - 1);
t = [0.005, 0.02, 0.08]; % seconds
figure(5);
hold on
sigma(sys_u,[],3);
for i = 1:length(t)
    s = tf('s');
    actuatorTF = 1/((t(i)*s) + 1);
    delta = simplify(actuatorTF - 1);
    sigma(delta);
end
title('Small Gain Theorem Robustness Test')
legend('1 + inv(L)', 'tau = 0.005', 'tau = 0.02', 'tau = 0.08')
% xlabel('Frequency (rad/s)')
% ylabel('Magnitude (dB)')
hold off