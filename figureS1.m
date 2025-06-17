%% Supplemental Figure 1 - Quantile-Quantile Plot
%
% Required Files:
% - face_dataset_paper.csv
% 
% Required Toolboxes:
% - Statistics and Machine Learning Toolbox
% 
% Code tested on MATLAB R2024b

clear; clc; clf; % clear all
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');  % removes table header warning
nrnd = 10; % digits for rounding from data import

%% Import and Grab Data
dtab = readtable('face_dataset_paper.csv'); % read d13C data from table 
Cind1 = strcmp(dtab.CO2Treatment,'Ambient'); % indices CO2 treatment, ambient
Xamb = round(dtab.d13C(Cind1),nrnd);  % d13C from ambient samples

%% Plotting
h = qqplot(Xamb); %quantile-qunatile plot
h(1).Marker = 'o';
h(1).MarkerEdgeColor = 'b';
h(1).LineWidth = 1;
h(2).Color = 'k';
h(2).LineWidth= 5;
h(3).Color = 'k';
h(3).LineWidth= 5;
axis([-3 3 -32 -26]);
xlabel('Standard Normal Quantiles');
ylabel('Quantiles of Ambient $\delta^{13}$C','Interpreter','latex');
title('');
set(gca,'Fontsize',24)