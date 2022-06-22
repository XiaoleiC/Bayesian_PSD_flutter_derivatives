close all
clear
clc
dt = 1/1024; % sampling interval
df = 1/dt; % sampling frequency
M = 8; % divide the signal into M pieces
%% structural properties (umit: SI)
B_deck = 0.9; % width of the bridge deck (unit: meter)
w_h = 2*pi*1.857; % vertical circular frequency
w_a = 2*pi*3.000; % torsional circular frequency
m = 11.7; % weight along longitudinal dimension (unit: kg/m)
I = 1.7; % inertial moment along longitudinal dimension (unit: kg*m^2/m)
zeta_h = 0.0020; % vertical damping ratio
zeta_a = 0.0015; % torsional damping ratio
rho_air = 1.225; % air density
%% wind field parameter
U = 6.77; % mean wind speed in the wind tunnel test
K = [B_deck*w_h/U; B_deck*w_a/U]; % reduced frequency
V_reduced = 2*pi./K; % reduced velocity
%% load the file of measured buffeting displacement respoonse (unit: meter)
temp = load('demo.mat');
y = temp.disp;
y(:,1) = filloutliers(y(:,1),'linear','grubbs'); % not necessary if the 
% original signal is in good condition.
y(:,2) = filloutliers(y(:,2),'linear','grubbs'); % not necessary if the 
% original signal is in good condition.
%% Theodorsen solutions (used for initial guess of identified flutter 
% derivatives)
F = @(k) 1 - 0.165 ./ (1 + (0.0455 ./ k).^2) - 0.335 ./ (1 + (0.3 ./ k).^2);
G = @(k) -0.165 * (0.0455 ./ k) ./ (1 + (0.0455 ./ k).^2) - ...
    0.335 * (0.3 ./ k) ./ (1 + (0.3 ./ k).^2);
H1 = @(K) -pi * F(K/2) ./ K;
H2 = @(K) -0.25 * pi ./ K .* (1 + F(K/2) + 4*G(K/2)./K);
H3 = @(K) -pi ./ (K.^2) .* (F(K/2) - G(K/2) .* K / 4);
H4 = @(K) pi / 4 .* (1 + 4*G(K/2)./K);

A1 = @(K) 0.25 * pi * F(K/2) ./ K;
A2 = @(K) -pi ./ (16*K) .* (1 - 4*G(K/2)./K - F(K/2));
A3 = @(K) pi./((2*K).^2) .* (1/32*K.^2 + F(K/2) - 0.25*G(K/2) .* K);
A4 = @(K) -0.25*pi .* (G(K/2)./K);
%% modified initial guess of identified flutter derivatives
A1_identified = rho_air*B_deck^3*w_h/(I)*A1(K(1));
A2_identified = rho_air*B_deck^4*w_a/(I)*A2(K(2));
A3_identified = rho_air*B_deck^4*w_a^2/(I)*A3(K(2));
A4_identified = rho_air*B_deck^3*w_h^2/(I)*A4(K(1));
H1_identified = rho_air*B_deck^2*w_h/(m)*H1(K(1));
H2_identified = rho_air*B_deck^3*w_a/(m)*H2(K(2));
H3_identified = rho_air*B_deck^3*w_a^2/(m)*H3(K(2));
H4_identified = rho_air*B_deck^2*w_h^2/(m)*H4(K(1));
%% PSD of measured buffeting displacement response
[psd_y_M,f,psd_plot] = fft_transfer(df,y,M,0); % y is the measured buffeting 
% signal; df is the sampling frequency; fft_transfer is a function used to
% get the PSD of y.
k1 = 405;
k2 = 510;
k3 = 675;
k4 = 790;
disp(['The PSD information in ',num2str(f(k1)),'-',num2str(f(k2)),...
    'Hz and ',num2str(f(k3)),'-',num2str(f(k4)),'Hz will be used.']);
%% for optimization: x0 is the initial guess; lb is the lower bound; ub is 
% the higher bound
magnify_up = [15,5,5,5,7,15,5,5];
magnify_bottom = [5,15,5,15,10,5,10,15];
lb = zeros(1,12);
ub = zeros(1,12);
x0 = [A1_identified,A2_identified,A3_identified,A4_identified,...
    H1_identified,H2_identified,H3_identified,H4_identified,-3.5,-3.5,-3.5,-3.5];
lb(1:8) = x0(1:8) - magnify_bottom.*abs(x0(1:8));
ub(1:8) = x0(1:8) + magnify_up.*abs(x0(1:8));
lb(9:end) = [-6,-6,-6,-6];
ub(9:end) = [-2,-2,-2,-2];
%% optimization: theta_opt is the identified result; 
% simulated annealing is used, time-consuming (about 20min) but robust; 
% MCMC can also be used since we have the likelihood function, much more 
% efficient but less robust. 
options = optimoptions('simulannealbnd','PlotFcns',...
    {@saplotbestx,@saplotbestf,@saplotx,@saplotf},'MaxFunctionEvaluations', 100000);
[theta_opt,fval,exitFlag,output] = simulannealbnd(@(theta) calculate_NLLF_noncross(k1, k2, k3, k4, M, w_h, w_a, zeta_h, zeta_a, theta, psd_y_M, f),x0,lb,ub,options);
%% Result display
A_1 = [A1_identified;theta_opt(1);abs((theta_opt(1) - A1_identified)/A1_identified)*100;lb(1);ub(1)];
A_2 = [A2_identified;theta_opt(2);abs((theta_opt(2) - A2_identified)/A2_identified)*100;lb(2);ub(2)];
A_3 = [A3_identified;theta_opt(3);abs((theta_opt(3) - A3_identified)/A3_identified)*100;lb(3);ub(3)];
A_4 = [A4_identified;theta_opt(4);abs((theta_opt(4) - A4_identified)/A4_identified)*100;lb(4);ub(4)];
H_1 = [H1_identified;theta_opt(5);abs((theta_opt(5) - H1_identified)/H1_identified)*100;lb(5);ub(5)];
H_2 = [H2_identified;theta_opt(6);abs((theta_opt(6) - H2_identified)/H2_identified)*100;lb(6);ub(6)];
H_3 = [H3_identified;theta_opt(7);abs((theta_opt(7) - H3_identified)/H3_identified)*100;lb(7);ub(7)];
H_4 = [H4_identified;theta_opt(8);abs((theta_opt(8) - H4_identified)/H4_identified)*100;lb(8);ub(8)];

Info = table(A_1,A_2,A_3,A_4,H_1,H_2,H_3,H_4,'RowNames',{'Theoretical values';'Identified values';'variance (%)';'lower band';'upper band'});
disp(Info);
%% verification (reconstruct the PSDs of buffeting displacement response 
% with identified theta_opt and measured signals)
f = f(1:1.2*k4);
f_temp = f(1:1.2*k4);
psd_measured = zeros(length(f),2);
psd_theoretical = zeros(length(f),2);
psd_identified = zeros(length(f),2);
psd_buffeting_plot = zeros(length(f),2);
psd_fL1 = 10^theta_opt(9);
psd_fM1 = 10^theta_opt(10);
psd_fL2 = 10^theta_opt(11);
psd_fM2 = 10^theta_opt(12);
for k3_temp = 1:size(psd_measured,1)
    if (k3_temp <= k2) && (k3_temp >= k1)
        H_temp1 = H_omega(w_h, w_a, zeta_h, zeta_a, A1_identified, A2_identified, A3_identified, A4_identified,...
            H1_identified, H2_identified, H3_identified, H4_identified, 2*pi*f(k3_temp));
        H_temp2 = H_omega(w_h, w_a, zeta_h, zeta_a, theta_opt(1), theta_opt(2), theta_opt(3), theta_opt(4),...
            theta_opt(5), theta_opt(6), theta_opt(7), theta_opt(8), 2*pi*f(k3_temp));
        psd_temp1 = H_temp1*diag([psd_fL1,psd_fM1])*H_temp1';
        psd_temp2 = H_temp2*diag([psd_fL1,psd_fM1])*H_temp2';
        psd_buffeting_plot(k3_temp,1) = log10(psd_fL1);
        psd_buffeting_plot(k3_temp,2) = log10(psd_fM1);
        psd_measured(k3_temp,1) = log10(real(psd_y_M{k3_temp}(1,1)));
        psd_measured(k3_temp,2) = log10(real(psd_y_M{k3_temp}(2,2)));
        psd_identified(k3_temp,1) = log10(real(psd_temp2(1,1)));
        psd_identified(k3_temp,2) = log10(real(psd_temp2(2,2)));
        psd_theoretical(k3_temp,1) = log10(real(psd_temp1(1,1)));
        psd_theoretical(k3_temp,2) = log10(real(psd_temp1(2,2)));
    elseif (k3_temp >= k3) && (k3_temp <=k4)
        H_temp1 = H_omega(w_h, w_a, zeta_h, zeta_a, A1_identified, A2_identified, A3_identified, A4_identified,...
            H1_identified, H2_identified, H3_identified, H4_identified, 2*pi*f(k3_temp));
        H_temp2 = H_omega(w_h, w_a, zeta_h, zeta_a, theta_opt(1), theta_opt(2), theta_opt(3), theta_opt(4),...
            theta_opt(5), theta_opt(6), theta_opt(7), theta_opt(8), 2*pi*f(k3_temp));
        psd_temp1 = H_temp1*diag([psd_fL2,psd_fM2])*H_temp1';
        psd_temp2 = H_temp2*diag([psd_fL2,psd_fM2])*H_temp2';
        psd_buffeting_plot(k3_temp,1) = log10(psd_fL2);
        psd_buffeting_plot(k3_temp,2) = log10(psd_fM2);
        psd_measured(k3_temp,1) = log10(real(psd_y_M{k3_temp}(1,1)));
        psd_measured(k3_temp,2) = log10(real(psd_y_M{k3_temp}(2,2)));
        psd_identified(k3_temp,1) = log10(real(psd_temp2(1,1)));
        psd_identified(k3_temp,2) = log10(real(psd_temp2(2,2)));
        psd_theoretical(k3_temp,1) = log10(real(psd_temp1(1,1)));
        psd_theoretical(k3_temp,2) = log10(real(psd_temp1(2,2)));
    else
        H_temp1 = H_omega(w_h, w_a, zeta_h, zeta_a, A1_identified, A2_identified, A3_identified, A4_identified,...
            H1_identified, H2_identified, H3_identified, H4_identified, 2*pi*f(k3_temp));
        H_temp2 = H_omega(w_h, w_a, zeta_h, zeta_a, theta_opt(1), theta_opt(2), theta_opt(3), theta_opt(4),...
            theta_opt(5), theta_opt(6), theta_opt(7), theta_opt(8), 2*pi*f(k3_temp));
        psd_temp1 = H_temp1*diag([0.5*(psd_fL1+psd_fL2),0.5*(psd_fM1+psd_fM2)])*H_temp1';
        psd_temp2 = H_temp2*diag([0.5*(psd_fL1+psd_fL2),0.5*(psd_fM1+psd_fM2)])*H_temp2';
        psd_buffeting_plot(k3_temp,1) = 0;
        psd_buffeting_plot(k3_temp,2) = 0;
        psd_measured(k3_temp,1) = log10(real(psd_y_M{k3_temp}(1,1)));
        psd_measured(k3_temp,2) = log10(real(psd_y_M{k3_temp}(2,2)));
        psd_identified(k3_temp,1) = 0;
        psd_identified(k3_temp,2) = 0;
        psd_theoretical(k3_temp,1) = 0;
        psd_theoretical(k3_temp,2) = 0;
    end
end
figure
hold on
plot(f(1:k4),psd_measured(1:k4,1),'linewidth',1.5,'linestyle','-');
scatter(f_temp(psd_identified(:,1)~=0),psd_identified(psd_identified(:,1)~=0,1),'r');
xlabel('Hz');
ylabel('PSD (log10)');
legend('Measured PSD','Reconstructed PSD');
title(['U',num2str(U),': PSD of vertical displacement']);
figure
hold on
plot(f(1:k4),psd_measured(1:k4,2),'linewidth',1.5,'linestyle','-');
scatter(f_temp(psd_identified(:,2)~=0),psd_identified(psd_identified(:,2)~=0,2),'r');
legend('Measured PSD','Reconstructed PSD');
xlabel('Hz');
ylabel('PSD (log10)');
title(['U',num2str(U),': PSD of torsional displacement']);
%% calculate lambda
reconstructed_PSD = [psd_identified(psd_identified(:,1)~=0,1);...
    psd_identified(psd_identified(:,2)~=0,2)];
measured_PSD = [psd_measured(psd_identified(:,1)~=0,1);psd_measured(psd_identified(:,2)~=0,2)]; 
lambda0 = calculate_lambda(reconstructed_PSD, measured_PSD);
disp(['lambda is ',num2str(lambda0),'.']);
