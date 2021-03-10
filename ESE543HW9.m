%% Question 1 a-c
clear all
clc

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

% Creating wiggle matrix
Aw = [ 0 Cp
0.*ones(2,1) Ap];
Bw = [ 0 Bp']';

% Checking if Aw and Bw is controllable
Kalman = [Bw Aw*Bw Aw^2*Bw];
rankKal = rank(Kalman);
% Rank of Kalman matrix is 3 which is same size as row of Aw so system is
% controllable

% Desired poles
dp1 = -5;
dp2 = -6;
dp3 = -7;
% Placing poles
Kp = place(Aw, Bw, [dp1, dp2, dp3]);
% Creating new wiggle matrix with desired poles
nAw = Aw-(Bw*Kp);
eigNAw = eig(nAw);
% Eigenvalues are the desired poles so pole placement was done correctly

% Creating linear quadratic regulator(LQR)
R = 1;
Q=0.*nAw;
Q(1,1)=10.;
[K,S,E]=lqr(nAw,Bw,Q,R);

% Creating Closed Loop System
% xdot = (Aw-Bw*K)x + Fr
% ycl = x(2)
F=[-1 0 0]';
Acl = nAw-Bw*K;
Bcl = F;
Ccl=[0 1 0];
Dcl = 0.*Ccl*F;

% Simulating system at constant 10 deg command
time = 0:0.1:10;
r = ones(1,length(time));
% Closed loop state space
SSCL = ss(Acl, Bcl, Ccl, Dcl);
[ycl,xcl] = lsim(SSCL,10*r,time);
% Plotting
figure(1)
hold on
plot(xcl,10*r)
plot(xcl,ycl)
grid on
title('Step Response of Angle of Attack at 10 degrees')
xlabel('Time (sec)')
ylabel('Y (states) and R (command)  Unit (degrees)')
legend('command','states')
hold off

% Common Controler Form
Ac = [0];
Bc1 = [1 0];
Bc2 = [-1];
Cc = -[K(1)];
Dc1 = -[K(2) K(3)];
Dc2 = [0];

% Closed Loop System
%State Space for Closed Loop System
Cp = [1 0; 0 0];
Dp = 0.*Cp*Bp;
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl1 = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl1 = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl1 = Bcl(:,1);
Ccl1 = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl1 =(Dp*Z*Dc2);
Dcl1 = Dcl(:,1);
syscl = ss(Acl1,Bcl1,Ccl1,Dcl1);
% Lu
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

% Nyquist for Lu
figure(2);
nyquist(sys_u)
axis([-2 2 -2 2]);

% Bode for Lu
figure(3);
margin(sys_u)
freqResp = allmargin(sys_u);

% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));
figure(4)
sigma(sys_u,[],2);
title('Singular Values for I + Lu')
% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));
figure(5)
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

% Plotting Sensitivity and Comp Sensitivity
[nCp,nAp] = size(Cp);
w = logspace(-1,3,500);
T  = freqresp(syscl,w);
S = 1 - T;
for jj = 1:nCp
    figure('Name','Co-Sensitivity (T)'),
    semilogx(w,20*log10(abs(squeeze(T(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(T(jj,1,:)))));
    legend([ 'T(' num2str(jj) ') max = ' num2str(v_min)],'Location','Best');
    %legend('T-Az LQR SF','T-Az LQG OF','T-Az LQR OFW','Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Co-Sensitivity (T)')
end

for jj = 1:nCp
    figure('Name','Sensitivity (S)'),
    semilogx(w,20*log10(abs(squeeze(S(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(S(jj,1,:)))));
    legend([num2str(jj) ' S max = ' num2str(v_min)],'Location','Best');
    %legend('T-Az LQR SF','T-Az LQG OF','T-Az LQR OFW','Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Sensitivity (S)')
end

%% Question 1 d-e
clear all
clc

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
[K,S,E]=lqr(nAw,Bw,Q,R);

% Creating Closed Loop System
% xdot = (Aw-Bw*K)x + Fr
% ycl = x(2)
F=[0 -1 0 0]';
Acl = nAw-Bw*K;
Bcl = F;
Ccl=[0 0 1 0];
Dcl = 0.*Ccl*F;

% Simulating system at constant 10 deg command
time = 0:0.1:10;
r = sin(w0*time);
% Closed loop state space
SSCL = ss(Acl, Bcl, Ccl, Dcl);
[ycl,xcl] = lsim(SSCL,r,time);
% Plotting
figure(1)
hold on
plot(xcl,r)
plot(xcl,ycl)
grid on
title('Autopilot Tracking Sinusoidal Command for Alpha')
xlabel('Time (sec)')
ylabel('Y (states) and R (command)  Unit (degrees)')
legend('command','states')
hold off

% Common Controler Form
Ac = [0 1; -w0^2 0];
Bc1 = [0 0; 1 0];
Bc2 = [0 -1]';
Cc = -[K(2) K(1)];
Dc1 = -[K(3) K(4)];
Dc2 = [0];

% Closed Loop System
%State Space for Closed Loop System
Cp = [1 0; 0 0];
Dp = 0.*Cp*Bp;
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl1 = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl1 = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl1 = Bcl(:,1);
Ccl1 = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl1 =(Dp*Z*Dc2);
Dcl1 = Dcl(:,1);
syscl = ss(Acl1,Bcl1,Ccl1,Dcl1);
% Lu
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

% Nyquist for Lu
figure(2);
nyquist(sys_u)
axis([-2 2 -2 2]);

% Bode for Lu
figure(3);
margin(sys_u)
freqResp = allmargin(sys_u);

% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));
figure(4)
sigma(sys_u,[],2);
title('Singular Values for I + Lu')
% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));
figure(5)
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

% Plotting Sensitivity and Comp Sensitivity
[nCp,nAp] = size(Cp);
w = logspace(-1,3,500);
T  = freqresp(syscl,w);
S = 1 - T;
for jj = 1:nCp
    figure('Name','Co-Sensitivity (T)'),
    semilogx(w,20*log10(abs(squeeze(T(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(T(jj,1,:)))));
    legend([ 'T(' num2str(jj) ') max = ' num2str(v_min)],'Location','Best');
    %legend('T-Az LQR SF','T-Az LQG OF','T-Az LQR OFW','Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Co-Sensitivity (T)')
end

for jj = 1:nCp
    figure('Name','Sensitivity (S)'),
    semilogx(w,20*log10(abs(squeeze(S(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(S(jj,1,:)))));
    legend([num2str(jj) ' S max = ' num2str(v_min)],'Location','Best');
    %legend('T-Az LQR SF','T-Az LQG OF','T-Az LQR OFW','Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Sensitivity (S)')
end

%% Question 2
clear all
clc

% Can't find magnetically suspended ball equations or model so using
% linearized suspended ball model found in HW#8
Ap = [0 1; 1 0];
Bp = [0 1]';
Cp = [1 0];
Dp = 0.*Cp*Bp;

% Creating wiggle system
Aw = [ 0 Cp
0.*ones(2,1) Ap];
Bw = [ 0 Bp']';

% Checking if Aw and Bw is controllable
Kalman = [Bw Aw*Bw Aw^2*Bw];
rankKal = rank(Kalman);
% Rank of Kalman matrix is 3 which is same size as row of Aw so system is
% controllable

% Desired poles
dp1 = -15;
dp2 = -20;
dp3 = -25;
% Placing poles
Kp = place(Aw, Bw, [dp1, dp2, dp3]);
% Creating new wiggle matrix with desired poles
nAw = Aw-(Bw*Kp);
eigNAw = eig(nAw);
% Eigenvalues are the desired poles so pole placement was done correctly

% Creating linear quadratic regulator(LQR)
R = 1;
Q=0.*nAw;
Q(1,1)=10.;
[K,S,E]=lqr(nAw,Bw,Q,R);

% Creating Closed Loop System
% xdot = (Aw-Bw*K)x + Fr
% ycl = x(2)
F=[-1 0 0]';
Acl = nAw-Bw*K;
Bcl = F;
Ccl=[0 1 0];
Dcl = 0.*Ccl*F;

% Simulating system at constant 10 deg command
time = 0:0.1:10;
r = ones(1,length(time));
% Closed loop state space
SSCL = ss(Acl, Bcl, Ccl, Dcl);
[ycl,xcl] = lsim(SSCL,10*r,time);
% Plotting
figure(1)
hold on
plot(xcl,10*r)
plot(xcl,ycl)
grid on
title('Step Response of Suspended Ball')
xlabel('Time (sec)')
ylabel('Y (states) and R (command)  Unit (meters)')
legend('command','states')
hold off

% Common Controler Form
Ac = [0];
Bc1 = [1 0];
Bc2 = [-1];
Cc = -[K(1)];
Dc1 = -[K(2) K(3)];
Dc2 = [0];

% Closed Loop System
%State Space for Closed Loop System
Cp = [1 0; 0 0];
Dp = 0.*Cp*Bp;
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl1 = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl1 = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl1 = Bcl(:,1);
Ccl1 = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl1 =(Dp*Z*Dc2);
Dcl1 = Dcl(:,1);
syscl = ss(Acl1,Bcl1,Ccl1,Dcl1);
% Lu
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

% Nyquist for Lu
figure(2);
nyquist(sys_u)
axis([-2 2 -2 2]);

% Bode for Lu
figure(3);
margin(sys_u)
freqResp = allmargin(sys_u);

% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));
figure(4)
sigma(sys_u,[],2);
title('Singular Values for I + Lu')
% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));
figure(5)
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

% Plotting Sensitivity and Comp Sensitivity
[nCp,nAp] = size(Cp);
w = logspace(-1,3,500);
T  = freqresp(syscl,w);
S = 1 - T;
for jj = 1:nCp
    figure('Name','Co-Sensitivity (T)'),
    semilogx(w,20*log10(abs(squeeze(T(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(T(jj,1,:)))));
    legend([ 'T(' num2str(jj) ') max = ' num2str(v_min)],'Location','Best');
    %legend('T-Az LQR SF','T-Az LQG OF','T-Az LQR OFW','Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Co-Sensitivity (T)')
end

for jj = 1:nCp
    figure('Name','Sensitivity (S)'),
    semilogx(w,20*log10(abs(squeeze(S(jj,1,:)))),'b','LineWidth',2);grid
    v_min = max(max(abs(squeeze(S(jj,1,:)))));
    legend([num2str(jj) ' S max = ' num2str(v_min)],'Location','Best');
    %legend('T-Az LQR SF','T-Az LQG OF','T-Az LQR OFW','Location','Best');
    xlabel('Frequency (rps)')
    ylabel('Mag (dB)')
    title('Sensitivity (S)')
end