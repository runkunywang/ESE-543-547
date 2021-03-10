%%Question 1

%Constants
ZaV = -1.3046; %(1/s) 
ZdV = -0.2142; %(1/s) 
Ma = 47.7109; %(1/s2) 
Md = -104.8346; %(1/s2)
V = 886.78; %(fps)
Za = ZaV*V;
Zd = ZdV*V;
%Controller Constants
Ka=-0.0015; 
Kq=-0.32; 
az=2.0; 
aq=6.0;

%1.1
%State Space for Plant
Ap = [ZaV 1; Ma 0];
Bp = [ZdV; Md];
Cp = [Za 0; 0 1];
Dp = [Zd; 0];
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
%checking eigenvalues
eAcl = eig(Acl);
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

%1.2
%Step response plot
figure(1);
step(syscl)
%63% rise time
[stepcl, tOut] = step(syscl);
stepRep = [stepcl tOut];
%Using data from stepRep, 63% rise time occurs at
rise63 = 0.3822; %Seconds
%Using data from stepRep, 95% settling time occurs at
setTime95 = 0.9750; %Seconds

%1.3
%Nyquist plot for Lu
figure(2);
nyquist(sys_u)
axis([-2 2 -2 2]);

%Nyquist plot for q in Ly with Az closed
%Nyquist plot for Az in Ly with Az closed
figure(3);
Azclose = sys_y*[1 0; 0 0];
nyquist(sys_y)
axis([-2 2 -2 2]);

%1.4 Lu Bode plot, LGCF, Phase crossover frequency
figure(4);
margin(sys_u)
[gm, pm, Wcg, Wcp] = margin(sys_u);
% LGCF = 4.53 rad/sec
% Phase Crossover Frequency = 32.9 rad/sec

%1.5 sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));
% minimum singular value is 0.9088
figure(5);
sigma(sys_u,[],2)

%1.6 sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));
% minimum singular value is 0.7389
figure(6);
sigma(sys_u,[],3)

%1.7 sigma(I+Ly)
sbar3 = sigma(sys_y,[],2);
minsbar3 = min(sbar3(1,:));
minsbar4 = min(sbar3(2,:));
% minimum singular value for Az is 60.8383
% minimum singular value for q is 0.0024
figure(7);
sigma(sys_y,[],2)

%1.8 Classical stability margins at plant input and singular value margins
gainMLu = 20*log10(gm); %-12.1 dB
phaseMLu = pm; %71.1 deg
RDu_nGM = 1/(1+minsbar1);
RDu_pGM = 1/(1-minsbar1);
RDu_Pha = 2*asin(minsbar1/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi ;
%1.9 Classical stability margins at plant output for q and Az
lyM = allmargin(sys_y);
figure(8);
bode(sys_y)

%1.10 Complementary Sensitivity T and Sensitivity S of Az
w = logspace(-3, 3, 500);
AzSS = sys_y*[1 0; 0 0];
T = freqresp(AzSS, w);
T = squeeze(T);
S = 1 - T;
Tproc1 = zeros(1,500);
Tproc2 = zeros(1,500);
Sproc1 = zeros(1,500);
Sproc2 = zeros(1,500);
for i = 1:500
    Tproc1(i) = T(1,1,i);
    Tproc2(i) = T(2,1,i);
    Sproc1(i) = S(1,1,i);
    Sproc2(i) = S(2,1,i);
end 
%Finding peak of T and S in dB
maxT1 = 20*log10(max(Tproc1));
maxT2 = 20*log10(max(Tproc2));
maxS1 = 20*log10(max(Sproc1));
maxS2 = 20*log10(max(Sproc2));
%Setting up w
w = logspace(-1,3,500);
figure(9)
hold on
semilogx(w,Tproc1)
set(gca, 'XScale', 'log');
title('Complementary Sensitivity (T)')
xlabel('Frequency')
ylabel('Complementary Sensitivity (Az/Azcmd)')
hold off

figure(10)
hold on
semilogx(w,Sproc1)
set(gca, 'XScale', 'log');
title('Sensitivity (S)')
xlabel('Frequency')
ylabel('Sensitivity (eAz/Azcmd)')
hold off



%% Question 2a

%Constants
ZaV = -1.3046; %(1/s) 
ZdV = -0.2142; %(1/s) 
Ma = 47.7109; %(1/s2) 
Md = -104.8346; %(1/s2)
V = 886.78; %(fps)
Za = ZaV*V;
Zd = ZdV*V;
%Controller Constants
Ka=-0.0015; 
Kq=-0.32; 
az=2.0; 
aq=6.0;
% Actuator Constants
%Natural Frequency
w_act = 150;
% Time constant
z_act = 0.6;

%2a.1 plant state space
Ap2 = [ZaV, 1, ZdV, 0; Ma, 0, Md, 0; 0, 0, 0, 1; 0, 0, -w_act*w_act, -2*z_act*w_act];
Bp2 = [0; 0; 0; w_act*w_act];
Cp2 = [Za, 0, Zd, 0; 0, 1, 0, 0];
Dp2 = [0; 0];
%State Space for Controller
Ac = [0 0; Kq*aq 0];
Bc1 = [-Ka*az 0; -Ka*Kq*aq -Kq*aq];
Bc2 = [Ka*az; Ka*Kq*aq];
Cc = [Kq 1];
Dc1 = [-Ka*Kq -Kq];
Dc2 = [Ka*Kq];
%State Space for Closed Loop System
Z = inv(eye(size(Dc1*Dp2))-Dc1*Dp2);  
Acl = [(Ap2+Bp2*Z*Dc1*Cp2) (Bp2*Z*Cc);(Bc1*(Cp2+Dp2*Z*Dc1*Cp2)) (Ac+Bc1*Dp2*Z*Cc)];

%checking eigenvalues
eAcl = eig(Acl);
Bcl = [Bp2*Z*Dc2;(Bc2+Bc1*Dp2*Z*Dc2)];
Bcl = Bcl(:,1);
Ccl = [(Cp2+Dp2*Z*Dc1*Cp2) (Dp2*Z*Cc)];
Dcl =(Dp2*Z*Dc2);
Dcl = Dcl(:,1);
syscl = ss(Acl,Bcl,Ccl,Dcl);
%SS model of loop gain Lu at the plant input
Ain= [Ap2 0.*Bp2*Cc; Bc1*Cp2 Ac];
Bin = [Bp2; Bc1*Dp2];
Cin= -[Dc1*Cp2 Cc];%change sign for loop gain
Din = -[Dc1*Dp2];
sys_u= ss(Ain,Bin,Cin,Din);
%SS model of loop gain Ly at the plant output
Aout= [Ap2 Bp2*Cc; 0.*Bc1*Cp2 Ac];
Bout = [Bp2*Dc1; Bc1];
Cout= -[Cp2 Dp2*Cc];%change sign for loop gain
Dout= -[Dp2*Dc1];
sys_y= ss(Aout,Bout,Cout,Dout);

%2a.2
%Step response plot
figure(1);
step(syscl)
%63% rise time
[stepcl, tOut] = step(syscl);
stepRep = [stepcl tOut];
%Using data from stepRep, 63% rise time occurs at
rise63 = 0.3809; %Seconds
%Using data from stepRep, 95% settling time occurs at
setTime95 = 0.9829; %Seconds

%2a.3
%Nyquist plot for Lu
figure(2);
nyquist(sys_u)
axis([-2 2 -2 2]);

%Nyquist plot for q in Ly with Az closed
%Nyquist plot for Az in Ly with Az closed
figure(3);
Azclose = sys_y*[1 0; 0 0];
nyquist(sys_y)
axis([-2 2 -2 2]);

%2a.4 Lu Bode plot, LGCF, Phase crossover frequency
figure(4);
margin(sys_u)
[gm, pm, Wcg, Wcp] = margin(sys_u);
% LGCF = 118 rad/sec
% Phase Crossover Frequency = 33.4 rad/sec

%2a.5 sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));
% minimum singular value is 0.6515
figure(5);
sigma(sys_u,[],2)

%2a.6 sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));
% minimum singular value is 0.7330
figure(6);
sigma(sys_u,[],3)

%2a.7 sigma(I+Ly)
sbar3 = sigma(sys_y,[],2);
minsbar3 = min(sbar3(1,:));
minsbar4 = min(sbar3(2,:));
% minimum singular value for Az is 55.49
% minimum singular value for q is 55.49
figure(7);
sigma(sys_y,[],2)

%2a.8 Classical stability margins at plant input and singular value margins
gainMLu = 20*log10(gm); %10.7 dB
phaseMLu = pm; %55.5 deg
RDu_nGM = 1/(1+minsbar1);
RDu_pGM = 1/(1-minsbar1);
RDu_Pha = 2*asin(minsbar1/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi ;
%2a.9 Classical stability margins at plant output for q and Az
lyM = allmargin(sys_y);
figure(8)
bode(sys_y)

%2a.10 Complementary Sensitivity T and Sensitivity S of Az
w = logspace(-3, 3, 500);
AzSS = sys_y*[1 0; 0 0];
T = freqresp(AzSS, w);
T = squeeze(T);
S = 1 - T;
Tproc1 = zeros(1,500);
Tproc2 = zeros(1,500);
Sproc1 = zeros(1,500);
Sproc2 = zeros(1,500);
for i = 1:500
    Tproc1(i) = T(1,1,i);
    Tproc2(i) = T(2,1,i);
    Sproc1(i) = T(1,1,i);
    Sproc2(i) = T(2,1,i);
end 
%Finding peak of T and S in dB
maxT1 = 20*log10(max(Tproc1));
maxT2 = 20*log10(max(Tproc2));
maxS1 = 20*log10(max(Sproc1));
maxS2 = 20*log10(max(Sproc2));
%Setting up w
w = logspace(-1,3,500);
figure(9);
hold on
semilogx(w,Tproc1)
set(gca, 'XScale', 'log');
title('Complementary Sensitivity (T)')
xlabel('Frequency')
ylabel('Complementary Sensitivity (Az/Azcmd)')
hold off

figure(10);
hold on
semilogx(w,Sproc1)
set(gca, 'XScale', 'log');
title('Sensitivity (S)')
xlabel('Frequency')
ylabel('Sensitivity (eAz/Azcmd)')
hold off

%% Question 2b
%Constants
ZaV = -1.3046; %(1/s) 
ZdV = -0.2142; %(1/s) 
Ma = 47.7109; %(1/s2) 
Md = -104.8346; %(1/s2)
V = 886.78; %(fps)
Za = ZaV*V;
Zd = ZdV*V;
%Controller Constants
Ka=-0.0015; 
Kq=-0.32; 
az=2.0; 
aq=6.0;
% Actuator Constants
%Natural Frequency
w_act = linspace(39,40,120);
% Time constant
z_act = 0.6;
% isstable evaluation matrix
stabM = zeros(1,120);
for i = 1:numel(stabM)
    %2a.1 plant state space
    Ap2 = [ZaV, 1, ZdV, 0; Ma, 0, Md, 0; 0, 0, 0, 1; 0, 0, -w_act(i)*w_act(i), -2*z_act*w_act(i)];
    Bp2 = [0; 0; 0; w_act(i)*w_act(i)];
    Cp2 = [Za, 0, Zd, 0; 0, 1, 0, 0];
    Dp2 = [0; 0];
    %State Space for Controller
    Ac = [0 0; Kq*aq 0];
    Bc1 = [-Ka*az 0; -Ka*Kq*aq -Kq*aq];
    Bc2 = [Ka*az; Ka*Kq*aq];
    Cc = [Kq 1];
    Dc1 = [-Ka*Kq -Kq];
    Dc2 = [Ka*Kq];
    %State Space for Closed Loop System
    Z = inv(eye(size(Dc1*Dp2))-Dc1*Dp2);  
    Acl = [(Ap2+Bp2*Z*Dc1*Cp2) (Bp2*Z*Cc);(Bc1*(Cp2+Dp2*Z*Dc1*Cp2)) (Ac+Bc1*Dp2*Z*Cc)];
    Bcl = [Bp2*Z*Dc2;(Bc2+Bc1*Dp2*Z*Dc2)];
    Bcl = Bcl(:,1);
    Ccl = [(Cp2+Dp2*Z*Dc1*Cp2) (Dp2*Z*Cc)];
    Dcl =(Dp2*Z*Dc2);
    Dcl = Dcl(:,1);
    syscl = ss(Acl,Bcl,Ccl,Dcl);
    stabM(i) = isstable(syscl);
end
finalM = [w_act; stabM]; 
%SS model of loop gain Lu at the plant input
Ain= [Ap2 0.*Bp2*Cc; Bc1*Cp2 Ac];
Bin = [Bp2; Bc1*Dp2];
Cin= -[Dc1*Cp2 Cc];%change sign for loop gain
Din = -[Dc1*Dp2];
sys_u= ss(Ain,Bin,Cin,Din);
%SS model of loop gain Ly at the plant output
Aout= [Ap2 Bp2*Cc; 0.*Bc1*Cp2 Ac];
Bout = [Bp2*Dc1; Bc1];
Cout= -[Cp2 Dp2*Cc];%change sign for loop gain
Dout= -[Dp2*Dc1];
sys_y= ss(Aout,Bout,Cout,Dout);

%2a.2
%Step response plot
figure(1);
step(syscl)
%63% rise time
[stepcl, tOut] = step(syscl);
stepRep = [stepcl tOut];
%Using data from stepRep, 63% rise time occurs at
rise63 = 0.3809; %Seconds
%Using data from stepRep, 95% settling time occurs at
setTime95 = 0.9829; %Seconds

%% Question 3

%Constants
K1 = 5;
K2 = -10;
K = [K1 0; 0 K2];
%transfer function
syms s;
h11 = 3.04/s;
h12 = -278.2/(s*(s+6.02)*(s+30.3));
h21 = 0.052/s;
h22 = -206.6/(s*(s+6.02)*(s+30.3));
H = [h11 h22; h21 h22];
% 3a multivariable Nyquist theorum
Dmat = [eye(2,2)-K*H];
detD = simplify(det(Dmat));
s = tf('s');
detDtf = (2500*s^4 + 52800*s^3 - 924145*s^2 - 12096428*s + 77165100)/(5*s^2*(10*s + 303)*(50*s + 301));
figure(1);
nyquist(detDtf)
axis([-10 10 -10 10]);

%transfer function
s = tf('s');
h11 = 3.04/s;
h12 = -278.2/(s*(s+6.02)*(s+30.3));
h21 = 0.052/s;
h22 = -206.6/(s*(s+6.02)*(s+30.3));
H = [h11 h22; h21 h22];
KH = K*H;
%3b singular values of return diff matrix and stab robust matrix
%singular values for return diff
sbar1 = sigma(KH,[],2);
figure(2)
sigma(KH,[],2);
%singular values for stab robust
sbar2 = sigma(KH,[],3);
figure(3)
sigma(KH,[],3);
%gain and phase margins
%creating matrix for each singular value
svrd1 = zeros(1,64);
svrd2 = zeros(1,64);
svsr1 = zeros(1,64);
svsr2 = zeros(1,64);
for i = 1:64
    svrd1(i) = sbar1(1,i);
    svrd2(i) = sbar1(2,i);
    svsr1(i) = sbar2(1,i);
    svsr2(i) = sbar2(2,i);
end
%creating gain margin matrices
NegGM1 = zeros(1,64);
NegGM2 = zeros(1,64);
NegGM3 = zeros(1,64);
NegGM4 = zeros(1,64);
PosGM1 = zeros(1,64);
PosGM2 = zeros(1,64);
PosGM3 = zeros(1,64);
PosGM4 = zeros(1,64);
for i = 1:64
    NegGM1(i) = 20*log10(1/(1+svrd1(i)));
    NegGM2(i) = 20*log10(1/(1+svrd2(i)));
    NegGM3(i) = 20*log10(1/(1+svsr1(i)));
    NegGM4(i) = 20*log10(1/(1+svsr2(i)));
    PosGM1(i) = 20*log10(1/(1-svrd1(i)));
    PosGM2(i) = 20*log10(1/(1-svrd2(i)));
    PosGM3(i) = 20*log10(1/(1-svsr1(i)));
    PosGM4(i) = 20*log10(1/(1-svsr2(i)));
end
%Calculating Gain Margins
mNGM1 = min(NegGM1);
mNGM2 = min(NegGM2);
mNGM3 = min(NegGM3);
mNGM4 = min(NegGM4);
mPGM1 = min(PosGM1);
mPGM2 = min(PosGM2);
mPGM3 = min(PosGM3);
mPGM4 = min(PosGM4);

%phase margin matrices
PM1 = zeros(1,64);
PM2 = zeros(1,64);
PM3 = zeros(1,64);
PM4 = zeros(1,64);
for i = 1:64
    PM1(i) = 180*(2*asin(svrd1(i)/2))/pi;
    PM2(i) = 180*(2*asin(svrd2(i)/2))/pi;
    PM3(i) = 180*(2*asin(svsr1(i)/2))/pi;
    PM4(i) = 180*(2*asin(svsr2(i)/2))/pi;
end
%Calculating Phase Margins
mPM1 = min(PM1);
mPM2 = min(PM2);
mPM3 = min(PM3);
mPM4 = min(PM4);

%preparing w
% [sbarT, w] = sigma(KH,[],2);
% w = w';
w = logspace(-1,3,64);
%Plotting gain margins
figure(4)
hold on
semilogx(w,NegGM1)
semilogx(w,PosGM1)
set(gca, 'XScale', 'log');
title('Gain Margin for Singular Values of Return Difference Matrix')
xlabel('Frequency (rad/s)')
ylabel('Gain Margin (dB)')
legend({'y = Negative Gain Margin','y = Positive Gain Margin'},'Location','southwest')
hold off

figure(5)
hold on
semilogx(w,NegGM2)
semilogx(w,PosGM2)
set(gca, 'XScale', 'log');
title('Gain Margin for Singular Values of Return Difference Matrix')
xlabel('Frequency (rad/s)')
ylabel('Gain Margin (dB)')
legend({'y = Negative Gain Margin','y = Positive Gain Margin'},'Location','southwest')
hold off

figure(6)
hold on
semilogx(w,NegGM3)
semilogx(w,PosGM3)
set(gca, 'XScale', 'log');
title('Gain Margin for Singular Values of Stability Robustness Matrix')
xlabel('Frequency (rad/s)')
ylabel('Gain Margin (dB)')
legend({'y = Negative Gain Margin','y = Positive Gain Margin'},'Location','southwest')
hold off

figure(7)
hold on
semilogx(w,NegGM4)
semilogx(w,PosGM4)
set(gca, 'XScale', 'log');
title('Gain Margin for Singular Values of Stability Robustness Matrix')
xlabel('Frequency (rad/s)')
ylabel('Gain Margin (dB)')
legend({'y = Negative Gain Margin','y = Positive Gain Margin'},'Location','southwest')
hold off

%Plotting phase margin
figure(8)
hold on
semilogx(w,PM1)
set(gca, 'XScale', 'log');
title('Phase Margin for Singular Values of Return Difference Matrix')
xlabel('Frequency (rad/s)')
ylabel('Degrees')
hold off

figure(9)
hold on
semilogx(w,PM2)
set(gca, 'XScale', 'log');
title('Phase Margin for Singular Values of Return Difference Matrix')
xlabel('Frequency (rad/s)')
ylabel('Degrees')
hold off

figure(10)
hold on
semilogx(w,PM3)
set(gca, 'XScale', 'log');
title('Phase Margin for Singular Values of Stability Robustness Matrix')
xlabel('Frequency (rad/s)')
ylabel('Degrees')
hold off

figure(11)
hold on
semilogx(w,PM4)
set(gca, 'XScale', 'log');
title('Phase Margin for Singular Values of Stability Robustness Matrix')
xlabel('Frequency (rad/s)')
ylabel('Degrees')
hold off

