%% Generates Bar Plot in Figure 2
%
% Required Files:
% - model_data.csv
% 
% Required Toolboxes:
% - N/A
% 
% Code tested on MATLAB R2024b

clear; clc;

dtab = readtable('model_data.csv'); % read model data from McCormack (2015)
ind = [1 7]; % Century and Landcarb row indices in table

% Rates (in yr^-1)
rts = [0.1 5]'; % slow, fast, intermediate turnover rates (yr^-1)
dcr = [1.46; 1]; % DCDM rates, fine root and transport root (yr^-1)
mrt = [1./3.38; 2./(3.77+4.15)]; % rates estimated from minirhizotron data

%% Total Eco C
for i = 1:2
dat = table2array(dtab(ind(i):(ind(i)+1),3)); % convert data to array
a = polyfit(rts,dat,1); % linear interpolation of McCormack estimates

tec_rts = a(2) + a(1)*rts; % TEC estimates given McCormack rates
tec_dcr = a(2) + a(1)*dcr; % TEC estimates given DCDM rates
tec_mrt = a(2) + a(1)*mrt; % TEC estimates given minirhizitron rates

% as percentage away from minirhizotron
tec_rts_per(:,i) = tec_rts./tec_mrt - 1; % percentage from McCormack
tec_dcr_per(:,i) = tec_dcr./tec_mrt - 1; % percentage for DCDM 
tec_mrt_per(:,i) = tec_mrt./tec_mrt - 1; % minirhizitron, should be zero
end


%% Soil Organic C
for i = 1:2
dat = table2array(dtab(ind(i):(ind(i)+1),7)); % convert data to array
a = polyfit(rts,dat,1); % linear interpolation of McCormack estimates

soc_rts = a(2) + a(1)*rts; % SOC estimates given McCormack rates
soc_dcr = a(2) + a(1)*dcr; % SOC estimates given DCDM rates
soc_mrt = a(2) + a(1)*mrt; % SOC estimates given minirhizitron rates

% as percentage away from minirhizotron
soc_rts_per(:,i) = soc_rts./soc_mrt - 1; % percentage from McCormack
soc_dcr_per(:,i) = soc_dcr./soc_mrt - 1; % percentage for DCDM 
soc_mrt_per(:,i) = soc_mrt./soc_mrt - 1; % minirhizitron, should be zero
end

%% NPP
for i = 1:2

dat = table2array(dtab(ind(i):(ind(i)+1),8)); % convert data to array
a = polyfit(rts,dat,1);  % linear interpolation of McCormack estimates

npp_rts = a(2) + a(1)*rts; % NPP estimates given McCormack rates
npp_dcr = a(2) + a(1)*dcr; % NPP estimates given DCDM rates
npp_mrt = a(2) + a(1)*mrt; % NPP estimates given minirhizitron rates

% as percentage from minirhizotron
npp_rts_per(:,i) = npp_rts./npp_mrt - 1; % percentage from McCormack
npp_dcr_per(:,i) = npp_dcr./npp_mrt - 1; % percentage for DCDM 
npp_mrt_per(:,i) = npp_mrt./npp_mrt - 1; % minirhizitron, should be zero
end


%% Plotting
xstrng = categorical(["CENTURY", "LANDCARB"]);
plts_grp1 = [ tec_dcr_per(1,:); soc_dcr_per(1,:); npp_dcr_per(1,:)]; 
plts_grp2 = [ tec_dcr_per(2,:); soc_dcr_per(2,:); npp_dcr_per(2,:)];

figure(1);
subplot(1,2,1) % Orders 1,2
b1 = bar(xstrng,100*plts_grp1); hold on;
set(gca,'Fontsize',28)
legend('Ecosystem C', 'Soil C','NPP')
ylim([-20 20])

subplot(1,2,2) % Orders 3,4,5
b2 = bar(xstrng,100*plts_grp2); hold on;
set(gca,'Fontsize',28)
legend('Ecosystem C', 'Soil C','NPP')
ylim([-20 20])

