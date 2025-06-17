%% Figure 3 - Pseudocolor Plot of DCDM Distribution and Mean
%
% Required Files:
% - face_dataset_paper.csv
% 
% Required Toolboxes:
% - Statistics and Machine Learning Toolbox
% - Financial Toolbox
% - Curve Fitting Toolbox
% 
% Code tested on MATLAB R2024b

clc; clf; % clear command and figure windows
clearvars -except dtab; % clear variables except table
opts3 = odeset('RelTol',1e-10,'AbsTol',1e-12); % ODE options
nrnd = 10; % digits for rounding from data import

%% Import and Initialize Data
dtab = readtable('face_datasetf.xlsx'); % read d13C data table 
tshft = 18; % time offset (in days) of initial data from t=0
t0 = min(dtab.DateSampled)-tshft; % set initial time t=0 to be cessation of fumigation
t = daysact(t0,dtab.DateSampled)/365; % set time (in years) for samples

% Indexing
Cind1 = strcmp(dtab.CO2Treatment,'Ambient'); % indices CO2 treatment, ambient
Cind2 = strcmp(dtab.CO2Treatment,'Elevated'); % indices CO2 treatment, elevated

% Fixed Parameters 
T = 4; % end time (in years)
ts = 0:.01:T; % time grid
xeq = mean(dtab.d13C(Cind1));  % equilibrium d13C, estimated from ambient data

Ts = {3}; % initialize cell arrays for time values
Xs = {3}; % initialize cell arrays for d13C values
for i = 1:3
    if i == 1, Oind = (dtab.Order == 1) | (dtab.Order == 2); end % Orders 1 and 2
    if i == 2, Oind = (dtab.Order == 3) | (dtab.Order == 4); end % Orders 2 and 3 
    if i == 3, Oind = (dtab.Order == 5); end % Order 5

    Ts{i} = t(Cind2 & Oind); % cell array for time values for each category
    Xs{i} = round(dtab.d13C(Cind2 & Oind),nrnd);  % cell array for elevated d13C values for each order category
end

Tv = [Ts{1};Ts{2};Ts{3}]; % combined time values
Xv = [Xs{1};Xs{2};Xs{3}]; % combined d13C values


%% Run Model Comparisons
% Assume:
% - xeq is known, fit from control (ambient) data above
% - x0, rates, and standard deviation are free parameters

% Define model parameters - taken from Table 1
% single exponential model
r_SEM = [0.63 0.48 0.48]; % rate for each order category
x0_SEM = -38.77; % estimated initial d13C

% demographic carbon distribution model
r_1 = 1.16; % r1 for DCDM 
r_2 = [1.46 0.98 1.01]; % turnover rate r2 for DCDM
x0_DCDM = -38; % estimated initial d13C
sd_DCDM = sqrt(2.71); % estimated initial standard deviation

% Define figure parameters
cmap_lines = [0.8 0.6 1; 0 1 0.5;0 0 1]; % Line colors
cmap = cool(3); % data point colors

x = -50:.1:-20; % d13C grid

for i = 1:3 % Iterate through: Orders 1,2; Orders 3,4; Orders 5
    v = [r_1, r_2(i), sd_DCDM]; % DCDM parameters
    P0 = normpdf(x,x0_DCDM,sd_DCDM); % initial distribution 
    [s,P] = ode45(@(t,P) Pfun(x,t,P,v,x0_DCDM,xeq),[0 T],P0,opts3); % run DCDM
    Pf = interp1(s,P,ts,'spline'); % interpolate to uniform time values

    pltmSPi = xeq + (x0_SEM-xeq)*exp(-r_SEM(i)*ts); % mean of SEM model

    A = r_2(i)/(r_2(i) - r_1); % coefficient
    pltmDCDM = xeq + (x0_DCDM - xeq)*( A*exp(-r_1*ts) + (1-A)*exp(-r_2(i)*ts) ); % mean of DCDM

    % Plotting by Order
    subplot(1,3,i)
    hold on;
    pcolor(ts,x,Pf'); shading interp;
    plot(ts,pltmSPi,'Color',cmap_lines(2,:),'LineWidth',5)
    plot(ts,pltmDCDM,'Color',cmap_lines(3,:),'LineWidth',5);
    plot(Ts{i},Xs{i},'o','Markersize',8,'Color',[0 0 0],'MarkerFaceColor',cmap(i,:))

    colormap hot;
    axis([0 4 -42 -25])
    set(gca,'Fontsize',30)
    clim([0 0.2])
end

drawnow;
disp('Plotting Finished')


%% Auxillary Functions
% Diff Eq. for DCDM 
function dPdt = Pfun(x,t,P,v,x0,xeq)
    x = reshape(x,[],1); % forces column vector    
    r1 = v(1); % pull r1 parameter from v
    r2 = v(2); % pull r2 parameter from v
    sd = v(3); % pull standard deviation from v
    
    C = xeq + (x0-xeq)*exp(-r1*t); 
    dPdt = -r2*P + r2*normpdf(x,C,sd); % differential equation

end
