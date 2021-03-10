M = 2;
m = 0.1;
l = 0.5;
g = 9.81;
Ap = [0 1 0 0; (g*(M+m)/(M*l)) 0 0 0; 0 0 0 1; (-m*g/M) 0 0 0];
Bp = [0; -1/(M*l); 0; 1/M];
Cp = eye(4);
Dp = 0.*Cp*Bp;
%% 3.a 
% Finding open loop poles
syms s
I = eye(4);
olc = simplify(det(s*I-Ap));
poles = simplify(solve(olc == 0, s));

%% 3.b 
% Controllability
Pc = [Bp Ap*Bp Ap^2*Bp Ap^3*Bp];
rankPc = rank(Pc);
%full rank can convert to controllable canonical form
% PlantSys = ss(Ap,Bp,Cp,Dp);
% cPlantSys = canon(PlantSys, 'companion');
% [cAp, cBp, cCp, cDp] = ssdata(cPlantSys);
% CCFAp = transpose(cAp);
% CCFBp = transpose(cCp);
% CCFCp = transpose(cBp);
% CCFDp = cDp;
CCFAp = [0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 20.601 0];
CCFBp = [0; 0; 0; 1];

%% 3.c
% Finding T
PcBar = [CCFBp CCFAp*CCFBp CCFAp^2*CCFBp CCFAp^3*CCFBp];
invPcBar = inv(PcBar);
T = Pc*invPcBar;

%% 3.d
% Finding gains
% Desired Eigenvalues l1 = -10, l2 = -10, l3-4 = -2+/-j3.464
syms s k
clCE = expand((s+10)*(s+10)*(s - (-2 + 3.464j))*(s - (-2 - 3.464j)));
ABK = CCFAp - CCFBp*k;
k1 = 999956/625;
k2 = 2249956/3125;
% k3 = solve(20601/1000 - k == -562489/15625, k);
k3 = 216.6003;
k4 = 24;
Kbar = [k1 k2 k3 k4];
Kx =  Kbar*inv(T);

%% 3.e
% ii
% Controller State Space
Ac = 0;
Bc1 = [0 0 0 0];
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

% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));

% singular value stability margins
RDu_nGM = 1/(1+minsbar1);
RDu_pGM = 1/(1-minsbar1);
RDu_Pha = 2*asin(minsbar1/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi; 

%% 3.e.iii
% eigenvalue check
% if max(real(eig(Acl))) > 0
%     disp('Closed-Loop System is Unstable');
%     disp([' Most unstable eig = ', num2str(max(real(eig(Acl))))]);
% return
% else
%     disp('Closed-Loop System is stable');
%     damp(Acl)
% end

%% 3.e.iii
% simulation
x0 = [0.1; 0; 0; 0];
x_IC = [0.1; 0; 0; 0; 0]; 
time = 0:0.05:3;
u = ones(size(time));
figure(3)
lsim(syscl,u,time,x_IC);

%% Question 4
Ap = [0 1; 0 -10];
Bp = [0; 1];
Cp = [40 0];
Dp = 0.*Cp*Bp;

% Controllability
Pc = [Bp Ap*Bp];
rankPc = rank(Pc);

%% Finding gains
syms s l1 l2 
clCE = simplify(expand((s - (-15 + 15j))*(s - (-15 - 15j))));
Kx = [450 20];
clE = expand((s+50)^2);
L = [l1; l2];
ALC = Ap + L*Cp;
detALC = det(s*eye(2)-ALC);
L = -[-9/4; -40];

%% Controller State Space
Ac = Ap - Bp*Kx - L*Cp;
Bc1 = L;
Bc2 = 0.*L;
Cc = -Kx;
Dc1 = 0;
Dc2 = 0;

% Closed Loop State Space
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl = Bcl(:,1);
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
Dcl = Dcl(:,1);
syscl = ss(Acl,Bcl,Ccl,Dcl);


%% Step Response 1
x01 = [0; 0];
x02 = [0; 0];
x_IC = [x01; x02]; 
time = 0:0.05:20;
u = ones(size(time));
figure(1)
lsim(syscl,u,time,x_IC);

%% Step Response 2
x01 = [0; 0];
x02 = [0.5; 0];
x_IC = [x01; x02]; 
time = 0:0.05:20;
u = ones(size(time));
figure(2)
lsim(syscl,u,time,x_IC);

%% Lu 
% Lu
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

% Nyquist for Lu
figure(3);
nyquist(sys_u)
axis([-2 2 -2 2]);

% Bode for Lu
figure(4);
margin(sys_u)
freqResp = allmargin(sys_u);

% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));

% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));

% singular value stability margins
RDu_nGM = 1/(1+minsbar1);
RDu_pGM = 1/(1-minsbar1);
RDu_Pha = 2*asin(minsbar1/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi; 

%%
M = 2;
m = 0.1;
l = 0.5;
g = 9.81;
Ap = [0 1 0 0; (g*(M+m)/(M*l)) 0 0 0; 0 0 0 1; (-m*g/M) 0 0 0];
Bp = [0; -1/(M*l); 0; 1/M];
Cp = eye(4);
Dp = 0.*Cp*Bp;

%% 3.a 
% Finding open loop poles
syms s
I = eye(4);
olc = simplify(det(s*I-Ap));
poles = simplify(solve(olc == 0, s));

%% 3.b 
% Controllability
Pc = [Bp Ap*Bp Ap^2*Bp Ap^3*Bp];
rankPc = rank(Pc);
%full rank can convert to controllable canonical form
% PlantSys = ss(Ap,Bp,Cp,Dp);
% cPlantSys = canon(PlantSys, 'companion');
% [cAp, cBp, cCp, cDp] = ssdata(cPlantSys);
% CCFAp = transpose(cAp);
% CCFBp = transpose(cCp);
% CCFCp = transpose(cBp);
% CCFDp = cDp;
CCFAp = [0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 20.601 0];
CCFBp = [0; 0; 0; 1];

%% 3.c
% Finding T
PcBar = [CCFBp CCFAp*CCFBp CCFAp^2*CCFBp CCFAp^3*CCFBp];
invPcBar = inv(PcBar);
T = Pc*invPcBar;

%% 3.d
% Finding gains
% Desired Eigenvalues l1 = -10, l2 = -10, l3-4 = -2+/-j3.464
syms s k
clCE = expand((s+10)*(s+10)*(s - (-2 + 3.464j))*(s - (-2 - 3.464j)));
ABK = CCFAp - CCFBp*k;
k1 = 999956/625;
k2 = 2249956/3125;
% k3 = solve(20601/1000 - k == -562489/15625, k);
k3 = 216.6003;
k4 = 24;
Kbar = [k1 k2 k3 k4];
Kx =  Kbar*inv(T);

%% 3.e
% ii
% Controller State Space
Ac = 0;
Bc1 = [0 0 0 0];
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

% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));

% singular value stability margins
RDu_nGM = 1/(1+minsbar1);
RDu_pGM = 1/(1-minsbar1);
RDu_Pha = 2*asin(minsbar1/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi; 

%% 3.e.iii
% simulation
x0 = [0.1; 0; 0; 0];
x_IC = [0.1; 0; 0; 0; 0]; 
time = 0:0.05:3;
u = ones(size(time));
figure(3)
lsim(syscl,u,time,x_IC);

%%
%% Question 4
Ap = [0 1; 0 -10];
Bp = [0; 1];
Cp = [40 0];
Dp = 0.*Cp*Bp;

% Controllability
Pc = [Bp Ap*Bp];
rankPc = rank(Pc);

%% Finding gains
syms s l1 l2 
clCE = simplify(expand((s - (-15 + 15j))*(s - (-15 - 15j))));
Kx = [450 20];
clE = expand((s+50)^2);
L = [l1; l2];
ALC = Ap + L*Cp;
detALC = det(s*eye(2)-ALC);
L = -[-9/4; -40];

%% Controller State Space
Ac = Ap - Bp*Kx - L*Cp;
Bc1 = L;
Bc2 = 0.*L;
Cc = -Kx;
Dc1 = 0;
Dc2 = 0;

% Closed Loop State Space
Z = inv(eye(size(Dc1*Dp))-Dc1*Dp);
Acl = [(Ap+Bp*Z*Dc1*Cp) (Bp*Z*Cc);(Bc1*(Cp+Dp*Z*Dc1*Cp)) (Ac+Bc1*Dp*Z*Cc)];
Bcl = [Bp*Z*Dc2;(Bc2+Bc1*Dp*Z*Dc2)];
Bcl = Bcl(:,1);
Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)];
Dcl =(Dp*Z*Dc2);
Dcl = Dcl(:,1);
syscl = ss(Acl,Bcl,Ccl,Dcl);


%% Step Response 1
x01 = [0; 0];
x02 = [0; 0];
x_IC = [x01; x02]; 
time = 0:0.05:20;
u = ones(size(time));
figure(1)
lsim(syscl,u,time,x_IC);

%% Step Response 2
x01 = [0; 0];
x02 = [0.5; 0];
x_IC = [x01; x02]; 
time = 0:0.05:20;
u = ones(size(time));
figure(2)
lsim(syscl,u,time,x_IC);

%% Lu 
% Lu
Ain= [Ap 0.*Bp*Cc; Bc1*Cp Ac];
Bin = [Bp; Bc1*Dp];
Cin= -[Dc1*Cp Cc];%change sign for loop gain
Din = -[Dc1*Dp];
sys_u= ss(Ain,Bin,Cin,Din);

% Nyquist for Lu
figure(3);
nyquist(sys_u)
axis([-2 2 -2 2]);

% Bode for Lu
figure(4);
margin(sys_u)
freqResp = allmargin(sys_u);

% Singular Values
% sigma(I+Lu)
sbar1 = sigma(sys_u,[],2);
minsbar1 = min(sigma(sys_u,[],2));

% sigma(I+inv(Lu))
sbar2 = sigma(sys_u,[],3);
minsbar2 = min(sigma(sys_u,[],3));

% singular value stability margins
RDu_nGM = 1/(1+minsbar1);
RDu_pGM = 1/(1-minsbar1);
RDu_Pha = 2*asin(minsbar1/2);
RDu_nGM_dB = 20*log10(RDu_nGM);
RDu_pGM_dB = 20*log10(RDu_pGM);
RDu_Pha_deg = 180*RDu_Pha/pi; 