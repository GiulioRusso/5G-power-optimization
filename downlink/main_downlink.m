%% DOWNLINK
% hp. MIMO multiuser system (l x l) with Uniform Power Conditions and Power
% Optimization with Water Filling and Dinkelbach algorithm
%%

% Start script
clear all
close all 
clc

% Number of users
K = 10;

% Number of BR antennas
N = 64;

% Number of attempts
Np = 1000;

% Area
l = 300;
a = l^2;
X_cell = [-l/2:1:l/2];
Y_cell = [-l/2:1:l/2];
% Cell position
Uxc = 0;
Uyc = 0;

% Users random distribution (K x Np)
Ux = round(l.*rand(K,Np) - l/2);
Uy = round(l.*rand(K,Np) - l/2);

% MAT. (K x Np) Distance Recever - Transmiter for each random user from
% the cell
D = zeros(K, Np);
for np=1:Np
    for k=1:K
        D(k, np) = sqrt((Ux(k, np) - Uxc)^2 + (Uy(k, np) - Uyc)^2);
    end
end

% Path Loss and Power Noise parameters
PLo = 10^(-0.1 * 84);
do = 35;
No_dBm = -140;
No = (1e-3) * 10^(0.1 * No_dBm);
F = 1;
eta = 3.5;

% MAT. (K x Np) Path Loss for each random user
PL_R = zeros(K, Np);
for np=1:Np
    for k=1:K
        PL_R(k, np) = sqrt(2 * PLo)./sqrt(1 + (D(k, np)./do).^(eta));
    end
end

% MAT. (N x K x Np) Fast Fading
H = (randn(N, K, Np) + sqrt(-1)*randn(N, K, Np))/sqrt(2);

% MAT. (N x K x Np) Slow Fading
for np=1:Np
    h(:,:,np)  = PL_R(:,np)' .* H(:,:,np);
end

% SCA. Band
B = 1e9;

% SCA. Noise Power
Pn = No * B * F;

% SCA. Performance
Efficiency = 1;
Performance = 1 / Efficiency;

% VEC. Users Power in dB and Naturals
P_dB_max = [-30:2:10];
P_max = 10.^(0.1*P_dB_max);



% MAT. (K x length(P_max)) Uniform Power values
P = zeros(K, length(P_max));

% MAT. (K x length(P_max) x Np) of SNR, Channel State Information and Capacity
% for Uniform Power Conditions
SNR = zeros(K, length(P_max), Np);
CSI = zeros(K, length(P_max), Np);
C = zeros(K, length(P_max), Np);
% VEC. (length(P_max) x Np) Sum Capacity
Ctot = zeros(length(P_max), Np);

% MAT. (K x length(P_max) x Np) of SINR, Channel State Information and Capacity
% for Uniform Power Conditions with Interference
SINR = zeros(K, length(P_max), Np);
CSI_I = zeros(K, length(P_max), Np);
C_I = zeros(K, length(P_max), Np);
% VEC. (length(P_max) x Np) Sum Capacity
Ctot_I = zeros(length(P_max), Np);



% MAT. (K x length(P_max) x Np) Optimal Power values from Water Filling and
% Water Level
P_opt = zeros(K, length(P_max), Np);
WaterLevel = zeros(length(P_max), Np);

% MAT. (K x length(P_max) x Np) of SNR, Channel State Information and Capacity
% for Water Filling
SNR_WF = zeros(K, length(P_max), Np);
CSI_WF = zeros(K, length(P_max), Np);
C_WF = zeros(K, length(P_max), Np);
% VEC. (length(P_max) x Np) Sum Capacity
Ctot_WF = zeros(length(P_max), Np);

% MAT. (K x length(P_max) x Np) of SINR, Channel State Information and Capacity
% for Water Filling with Interference
SINR_WFI = zeros(K, length(P_max), Np);
CSI_WFI = zeros(K, length(P_max), Np);
C_WFI = zeros(K, length(P_max), Np);
% VEC. (length(P_max) x Np) Sum Capacity
Ctot_WFI = zeros(length(P_max), Np);



% VEC. (K x 1) Dissipated Power
P_c = 1;

% VEC. (length(P_max) x Np) Energy Efficiency
EE = zeros(length(P_max), Np);

% VEC. (K x length(P_max) x Np) Optimal Power values from Dinkelbach with Water Filling
P_opt_EE = zeros(K, length(P_max), Np);
% VEC. (length(P_max) x Np) Optimal Sum Capacity for EE
Ctot_EE = zeros(length(P_max), Np);
% VEC. (length(P_max) x Np) Energy Efficiency
EE_opt = zeros(length(P_max), Np);

% MAT. (K x length(P_max) x Np) of SINR, Channel State Information and Capacity
% for Dinkelbach with Interference
SINR_EE_opt_I = zeros(K, length(P_max), Np);
CSI_EE_opt_I = zeros(K, length(P_max), Np);
C_EE_opt_I = zeros(K, length(P_max), Np);
% VEC. (length(P_max) x Np) Sum Capacity
Ctot_EE_opt_I = zeros(length(P_max), Np);
% VEC. (length(P_max) x Np) Energy Efficiency
EE_opt_I = zeros(length(P_max), Np);

% VEC. (length(P_max) x Np) Energy Efficiency
EE_WF = zeros(length(P_max), Np);



% Sum Capacities and Energy Efficiency Calculation
for np=1:Np
    for i=1:length(P_max)
        
        % Uniform power conditions
        P(:,i) = (P_max(i) / K) * ones(K,1);
        % Rate
        [Ctot(i,np), C(:,i,np), SNR(:,i,np), CSI(:,i,np)] = SumCapacityCalc(h(:,:,np), Pn, P(:,i), B, false);
        [Ctot_I(i,np), C_I(:,i,np), SINR(:,i,np), CSI_I(:,i,np)] = SumCapacityCalc(h(:,:,np), Pn, P(:,i), B, true);
        % Energy Efficiency
        [EE(i,np)] = EnergyEfficiencyCalc(Ctot_I(i,np), Performance, P_max(i), P_c);
        
        % Power Optimization
        [P_opt(:,i,np), WaterLevel(i,np)] = WaterFilling(CSI(:,i,np), P_max(i));
        % Rate
        [Ctot_WF(i,np), C_WF(:,i,np), SNR_WF(:,i,np), CSI_WF(:,i,np)] = SumCapacityCalc(h(:,:,np), Pn, P_opt(:,i,np), B, false);
        [Ctot_WFI(i,np), C_WFI(:,i,np), SINR_WFI(:,i,np), CSI_WFI(:,i,np)] = SumCapacityCalc(h(:,:,np), Pn, P_opt(:,i,np), B, true);
        % Energy Efficiency
        [EE_WF(i,np)] = EnergyEfficiencyCalc(Ctot_I(i,np), Performance, sum(P_opt(:,i,np)), P_c);
        [EE_opt(i,np), P_opt_EE(:,i,np), Ctot_EE(i,np)] = Dinkelbach(B, CSI(:,i,np), Performance, P_c, Ctot(i,np), P_max(i), h(:,:,np), Pn, i, EE_WF(:,np), EE_opt(:,np), P_opt_EE(:,:,np));
        [Ctot_EE_opt_I(i,np), C_EE_opt_I(:,i,np), SINR_EE_opt_I(:,i,np), CSI_EE_opt_I(:,i,np)] = SumCapacityCalc(h(:,:,np), Pn, P_opt_EE(:,i,np), B, true);
        [EE_opt_I(i,np)] = EnergyEfficiencyCalc(Ctot_EE_opt_I(i,np), Performance, sum(P_opt_EE(:,i,np)), P_c);
        
    end
end

% Average value on all attempts for all Power values
Ctot_I_av = mean(Ctot_I, 2) / B; % (length(P_max) x Np) -> (length(P_max) x 1)
Ctot_WF_av = mean(Ctot_WF, 2) / B; % (length(P_max) x Np) -> (length(P_max) x 1)
Ctot_WFI_av = mean(Ctot_WFI, 2) / B; % (length(P_max) x Np) -> (length(P_max) x 1)
EE_av = mean(EE, 2) / 1e6; % (length(P_max) x Np) -> (length(P_max) x 1)
EE_opt_av = mean(EE_opt, 2) / 1e6; % (length(P_max) x Np) -> (length(P_max) x 1)
EE_opt_I_av = mean(EE_opt_I, 2) / 1e6; % (length(P_max) x Np) -> (length(P_max) x 1)

% Plot
sgtitle('DOWNLINK', 'FontName', 'Menlo');

% Random users plot
subplot(2, 4, [1 5]);
pbaspect([1 1 1]);
title('Random users distribution with the BR Station for one attempt');
hold on;
plot(Ux(:, Np), Uy(:, Np), 'x', 'linewidth', 1, 'markersize', 5);
axis ([-l/2 l/2 -l/2 l/2]);
xlabel ('m');
ylabel ('m');
% Cell plot
plot(Uxc, Uyc, '*', 'linewidth', 1, 'markersize', 5);
grid on;

% Sum Capacities plot
subplot(2, 4, [2 6]);
pbaspect([1 1 1]);
title(['Sum Rates average with antennas = ', num2str(N)]);
hold on;
% with WF and no interference
semilogy(P_dB_max, Ctot_WF_av, '-og', 'linewidth', 2, 'markersize', 3);
% with WF and interference
semilogy(P_dB_max, Ctot_WFI_av, '-or', 'linewidth', 2, 'markersize', 3);
% with Uniform Power Conditions and interference
semilogy(P_dB_max, Ctot_I_av, '-ob', 'linewidth', 2, 'markersize', 3);
xlabel ('Power dB [dB W]');
ylabel ('Sum Rate [bit / s / Hz]');
legend('No Interference and Water Filling','Interference and Water Filling', 'Interference and Uniform Power Conditions', 'Location', 'northwest');
grid on;

% Optimal Powers plot
subplot(2, 4, 3);
pbaspect([2 1 1]);
title(['Power distribution by Water Filling for one attempt with total Power = ', num2str(P_max(i))]);
hold on;
bar(P_opt(:,i,Np), 'r');
yline(WaterLevel(i,Np), '--');
xlim([0 K+1]);
xlabel ('Users');
ylabel ('1 / \nu');
grid on;
% CSI plot
subplot(2, 4, 7);
pbaspect([2 1 1]);
title('Inverse Channel State Information for one attempt');
hold on;
bar(1./CSI_I(:,i,Np));
xlim([0 K+1]);
xlabel ('Users');
ylabel ('1 / CSI');
grid on;

% Energy Efficiencies plot
subplot(2, 4, [4 8]);
pbaspect([1 1 1]);
title(['Energy Efficiencies average with antennas = ', num2str(N)]);
hold on;
% with Dinkelbach and no interference
semilogy(P_dB_max, EE_opt_av, '-o', 'linewidth', 2, 'markersize', 3);
% with Dinkelbach and interference
semilogy(P_dB_max, EE_opt_I_av, '-o', 'linewidth', 2, 'markersize', 3);
% with Uniform Power Conditions and interference
semilogy(P_dB_max, EE_av, '-o', 'linewidth', 2, 'markersize', 3);
xlabel ('Power dB [dB W]');
ylabel ('Energy Efficiency [Megabit / Joule]');
legend('No Interference and Dinkelbach', 'Interference and Dinkelbach', 'Interference and Uniform Power Conditions', 'Location', 'northwest');
grid on;

% End Script