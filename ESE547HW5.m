clear all
clc
close all

%% Other parameters
disp('****************************** Program Start ****************');
plot_file_name = 'Homework_5_OutputFeedback.ppt';
save_plots = 0; % Flag to bypass saving plots
w = logspace(-3,4,1000);
t = linspace(0,1.5,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);
rtd = 180/pi;

grav = 32.174;

m2ft = 3.2808;    % meters to feet conversion
ft2m = 1/3.2808;  % feet to meters conversion

%% Define Plant Matrices & Wiggle Matrix, Constants %%


% a - alpha, d - delta, w - omega, z - zeta
Za_V = -1.21;
Ma = 44.2506;
Zd_V = -.1987;
Md = -97.3213;
V = 886.78; % (fps)
Za = V*Za_V;
Zd = V*Zd_V;
w_act = 2*pi*11; % actuator natural frequency rps
z_act = 0.707; % actuator damping
%% Starting output feedback design
% State space A and B matrices for plant with acuator
Ap = [Za_V, Za, 0, Zd; 
      Ma/Za, 0, (Md-(Ma*Zd/Za)), 0; 
      0, 0, 0, 1; 
      0, 0, w_act*w_act, -2*w_act*z_act];

Bp = [0 0 0 w_act*w_act]';

Cp = eye(4);

Dp = 0.*Cp*Bp;

% Wiggle matrices for A and B
Aw = [0 1 0 0 0; zeros(4,1) Ap];

Bw = [0; Bp];

Cw = eye(5);

Dw = 0.*Cw*Bw;

ss_w = ss(Aw,Bw,Cw,Dw);

% Matrices to store result
riseTime = [];
stableMat = [];

% LQR
% penalty values
qq = logspace(-3, 2.5, 500);
Q=0.*Aw;
R = 1;

figure('Name','Acceleration Time History')
hold on
for i = 1:length(qq)
    Q(1,1)= qq(169);
%     Q(1,1) = 0.2448;
    [Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);
    
    % Controller implementation
    Ac_act = [0.];
    Bc1_act = [0. 1. 0. 0. 0.];
    Bc2_act = -1;
    Cc_act = [-Kx_lqr(1)];
    Dc1_act = [0 -Kx_lqr(2:5)];
    Dc2_act = 0.;

    Z = inv(eye(size(Dc1_act*Dw))-Dc1_act*Dw);
    Acl_act = [ (Aw+Bw*Z*Dc1_act*Cw) (Bw*Z*Cc_act);
    (Bc1_act*(Cw+Dw*Z*Dc1_act*Cw)) (Ac_act+Bc1_act*Dw*Z*Cc_act)];
    Bcl_act = [ Bw*Z*Dc2_act;
    (Bc2_act+Bc1_act*Dw*Z*Dc2_act)];
    Ccl_act = [(Cw+Dw*Z*Dc1_act*Cw) (Dw*Z*Cc_act)];
    Dcl_act =(Dw*Z*Dc2_act);
    sys_cl_act = ss(Acl_act,Bcl_act,Ccl_act,Dcl_act);
    
    % Output Feedback step response
    y = step(sys_cl_act,t);
    az = y(:,2); %  acceleration (fps2)
    aze = abs(ones(size(az))-az);  % error for az
    taur = 0.; taus= 0.; % rise time and settling time
    fv = aze(numel(aze)); % final value of the error
    e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taur = crosst(e_n1,t); % rise time 
    e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze));
    e_n1 = abs(e_n) + e_n;
    taus = crosst(e_n1,t); % settling time
    
    riseTime = [riseTime, [i; qq(i); taur; max(az); min(az)]];
    stableMat = [stableMat, [i; qq(i); isstable(sys_cl_act)]];
    if i == 169
        plot(t,az);grid
        Kpick = Kx_lqr;
    end
end

hold off

%% 

F = Aw-Bw*Kpick;
[FeigVec,FeigVal] = eig(F)

Cfeedback = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0];
Xr = [FeigVec(:,2:3), FeigVec(:,5)];
% Xr = FeigVec(:,3:5);
Ky = Kpick*Xr*inv(Cfeedback*Xr);

Ac_act2 = [0.];
Bc1_act2 = [0. 1. 0. 0. 0.];
Bc2_act2 = -1;
Cc_act2 = [-Ky(1)];
Dc1_act2 = [0. -Ky(2:3) 0 0];
Dc2_act2 = 0.;

TestMat = eye(6);
TestMat(1,1) = -1;

Z2 = inv(eye(size(Dc1_act2*Dw))-Dc1_act2*Dw);
Acl_act2 = [ (Aw+Bw*Z2*Dc1_act2*Cw) (Bw*Z2*Cc_act2);
            (Bc1_act2*(Cw+Dw*Z2*Dc1_act2*Cw)) (Ac_act2+Bc1_act2*Dw*Z2*Cc_act2)];
Bcl_act2 = [ Bw*Z2*Dc2_act2;
           (Bc2_act2+Bc1_act2*Dw*Z2*Dc2_act2)];
Ccl_act2 = [(Cw+Dw*Z2*Dc1_act2*Cw) (Dw*Z2*Cc_act2)];
Dcl_act2 =(Dw*Z2*Dc2_act2);
sys_cl_act2 = ss(Acl_act2,Bcl_act2,Ccl_act2,Dcl_act2);

% SS model of loop gain at the plant input Lu
A_Lu_act2 = [ Aw 0.*Bw*Cc_act2;  Bc1_act2*Cw Ac_act2];
B_Lu_act2 = [ Bw; Bc1_act2*Dw];
C_Lu_act2 = -[ Dc1_act2*Cw Cc_act2];%change sign for loop gain
D_Lu_act2 = -[ Dc1_act2*Dw];
sys_Lu_act2 = ss(A_Lu_act2,B_Lu_act2,C_Lu_act2,D_Lu_act2);
% OutputFeedbackResponse = allmargin(sys_Lu_act2) 
magdb_act2 = 20*log10(abs(squeeze(freqresp(sys_Lu_act2,w))));
wc2 = crosst(magdb_act2,w); % LGCF, assumes Lu is a scalar
sr_act2 = sigma(sys_Lu_act2,w,3);
sru_min_act2 = min(abs(sr_act2));
rd_act2 = sigma(sys_Lu_act2,w,2);
rdu_min_act2 = min(abs(rd_act2));
Lu_act2 = freqresp(sys_Lu_act2,w);
     
neg_gm =  min([ (1/(1+rdu_min_act2)) (1-sru_min_act2)]); % in dB
pos_gm =  max([ (1/(1-rdu_min_act2)) (1+sru_min_act2)]); % in dB
neg_gmdB = 20*log10( neg_gm ); % in dB
pos_gmdB = 20*log10( pos_gm ); % in dB
pm = 180*(max([2*asin(rdu_min_act2/2) 2*asin(sru_min_act2/2)]))/pi;% in deg
disp('Singular value margins')
disp(['Min Singular value I+Lu =    ' num2str(rdu_min_act2)])
disp(['Min Singular value I+invLu = ' num2str(sru_min_act2)])
disp(['Singular value gain margins = [' ...
         num2str(neg_gmdB) ' dB,' num2str(pos_gmdB) ' dB ]' ])
disp(['Singular value phase margins = [ +/-' ...
         num2str(pm)  ' deg ]' ])

% Analysis at plant output
T_act2  = freqresp(sys_cl_act2,w); % Complementary Sensitivity
S_act2 = 1 - T_act2; % Sensitivity
T_Az_act2 = 20*log10(abs(squeeze(T_act2(1,1,:))));
S_Az_act2 = 20*log10(abs(squeeze(S_act2(1,1,:))));
Tmax_act2 = max(T_Az_act2); % Inf Norm of T in dB
Smax_act2 = max(S_Az_act2); % Inf Norm of S in dB


Ain2 = A_Lu_act2;
Bin2 = B_Lu_act2;
Cin2 = C_Lu_act2;
Din2 = D_Lu_act2;


y2 = step(sys_cl_act,t);
az2 = y2(:,2); %  acceleration (fps2)
aze2 = abs(ones(size(az2))-az2);  % error for az
taur2 = 0.; taus2= 0.; % rise time and settling time
fv2 = aze2(numel(aze2)); % final value of the error
e_n2 = aze2 - fv2*ones(size(aze2)) - 0.36*ones(size(aze2));
e_n12 = abs(e_n2) + e_n2;
taur2 = crosst(e_n12,t); % rise time 
e_n2 = aze2 - fv2*ones(size(aze2)) - 0.05*ones(size(aze2));
e_n12 = abs(e_n2) + e_n2;
taus2 = crosst(e_n12,t); % settling time

figure('Name','Acceleration Time History with SPC output feedback control')
plot(t,az2,'r','LineWidth',2);grid
hold on
legend(['SPC output feedback response ' '63% Tr = ' num2str(taur2) ' 95% Ts = ' num2str(taus2)])
hold off
F2 = Aw-(Bw*[Ky 0 0]);
[FeigVec2,FeigVal2] = eig(F2);
FeigVal2(1,1) = -FeigVal2(1,1);
disp(diag(FeigVal2));
disp(FeigVec2);

% Nyquist plot
figure('Name','Nyquist Plot at Plant Input')
[Nyre, Nyim] = nyquist(sys_cl_act,w);
plotRe = squeeze(Nyre(2,:,:));
plotIm = squeeze(Nyim(2,:,:));
negRe = plotRe;
negIm = -plotIm;
        
hold on
plot(plotRe,plotIm,'b', 'LineWidth',2);grid
plot(negRe,negIm,'b', 'LineWidth',2)
plot(xx1,yy1,'r')
title('Nyquist Plot at Plant Input')
xlabel('Re(L)')
ylabel('Im(L)')
hold off

% Bode Plot
figure('Name', 'Bode Plot')
margin(sys_Lu_act2);grid

% I + lu Min Sing Val
figure('Name','Return Difference at Plant Input'),
semilogx(w,20*log10(abs(rd_act2)),'b','LineWidth',2); grid
legend(['SPC min(I+Lu) = ' num2str(rdu_min_act2)], ...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Return Difference at Plant Input')

figure('Name','Stability Robustness at Plant Input'),
semilogx(w,20*log10(abs(sr_act2)),'b', 'LineWidth',2);
legend(['SPC min(I+invLu) = ' num2str(sru_min_act2)], ...
    'Location','Best');
xlabel('Frequency (rps)')
ylabel('Mag dB')
title('Stability Robustness at Plant Input')

figure('Name','Comp Sens T');
semilogx(w,T_Az_act2,'g', 'LineWidth',2);grid
legend(['SPC ||T||inf = ' num2str(Tmax_act2) ' (dB)'],'Location','Best');
title('Comp Sens T');
xlabel('Freq (rps)');ylabel('Mag (dB)');

figure('Name','Sens S');
semilogx(w,S_Az_act2,'g', 'LineWidth',2);grid
legend(['SPC ||S||inf = ' num2str(Smax_act2) ' (dB)'],'Location','Best');
title('Sens S');
xlabel('Freq (rps)');ylabel('Mag (dB)');

% Compute the noise-to-control TF RSLQR
Z = inv(eye(size(Dc1_act2*Dw))-Dc1_act2*Dw);
Bv2 = [       Bw*Z*Dc1_act2;
     (Bc1_act2+Bc1_act2*Dw*Z*Dc1_act2)];
Cv2  = [ Z*Dc1_act2*Cw Z*Cc_act2];
Cvv2 = [ Cv2 ; Cv2*Acl_act2 ];
Dv2 = Z*Dc1_act2;
Dvv2 = [ Dv2; Cv2*Bv2];
sys_noise = ss(Acl_act2,Bv2,Cvv2,Dvv2);
v_2_u2  = freqresp(sys_noise,w); % Noise to control freq response
dele_Az2    = 20*log10(abs(squeeze(v_2_u2(1,1,:))));
dele_q2     = 20*log10(abs(squeeze(v_2_u2(1,3,:))));
deledot_Az2 = 20*log10(abs(squeeze(v_2_u2(2,1,:))));
deledot_q2  = 20*log10(abs(squeeze(v_2_u2(2,3,:))));

figure('Name','Noise-2-Control');
semilogx(w,dele_Az2,'k',w,dele_q2,'b',...
         'LineWidth',2);
legend('Az SPC','q SPC',...
       'Location','Best');
title(['Noise to Control ']);
xlabel('Freq (rps)');ylabel('Mag (dB)');
grid;

figure('Name','Noise-2-Control Rate');
semilogx(w,deledot_Az2,'k',w,deledot_q2,'b',...
         'LineWidth',2);
legend('Az SPC','q SPC',...
       'Location','Best');
title(['Noise to Control Rate']);
xlabel('Freq (rps)');ylabel('Mag (dB)');
grid;



