%% Supplemental Figure 2 - Variance Plots
%
% Required Files:
% - face_dataset_paper.csv
% 
% Required Toolboxes:
% - Optimization Toolbox
% - Financial Toolbox
% 
% Code tested on MATLAB R2024b

clc; clf; % clear command and figure windows
clearvars -except dtab; % clear variables except table
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');  % removes table header warning

opts1=  optimset('TolFun',1e-10,'TolX',1e-10,'display','off'); % least-squares display options
nrnd = 10; % digits for rounding from data import

%% Import and Initialize Data
dtab = readtable('face_dataset_paper.csv'); % read d13C data table 
tshft = 18; % time offset (in days) of data, largely negligible for fits
t0 = min(dtab.DateSampled)-tshft; % set inital time-sample as t = 0;
t = daysact(t0,dtab.DateSampled)/365; % set time (in years) for samples

% Indexing
Cind1 = strcmp(dtab.CO2Treatment,'Ambient'); % indices CO2 treatment, ambient
Cind2 = strcmp(dtab.CO2Treatment,'Elevated'); % indices CO2 treatment, elevated

% Fixed Parameters 
T = 4.1; % end time (in days)
ts = 0:0.01:T; 
xeq = mean(dtab.d13C(Cind1));  % mean equilibrium d13C estimated from ambient samples

%% Total Variance Calculuation
% Empirical Variance
Tv = t(Cind2);
Xv = round(dtab.d13C(Cind2),nrnd);
uts = unique(Tv); % find unique time points
vXv = zeros(1,length(uts)); % initialize variance array
cnt_ptsV = zeros(1,length(uts));
mu4V = zeros(1,length(uts));
for k = 1:length(uts)
    cnt_ptsV = length(Xv(Tv == uts(k)));
    vXv(k) = var( Xv(Tv == uts(k) ) ); % compute variance each time point
    mu4V(k) = mean( ( Xv(Tv == uts(k)) - mean( Xv(Tv == uts(k)) )  ).^4 );
end

seV = sqrt(  (mu4V - (cnt_ptsV-3).*(vXv.^2) ./(cnt_ptsV-1) )./(cnt_ptsV)  );

disp('Variance for Data (Total) Finished')

% DCDM Variance
V = VarFun(1.1674, 1.1342, -37.9896,xeq,1.6478, ts);
% Note: Values computed across total data instead of separated by order, 
% so these values will not exactly match the table values. See figure1.m

% f(t) Variance for Gaussian models
lbv = [-Inf -Inf -Inf];
ubv = [Inf Inf Inf];
f0 = [-2 10 0.5]; 
fcom = @(x) x(1).*uts.^2 + x(2).*uts + x(3) - vXv';
pvarcom = lsqnonlin(fcom,f0,lbv,ubv,opts1); 

% Generate variance at time points
fvar = @(x,t) x(1).*t.^2 + x(2).*t + x(3);
varGauss = fvar(pvarcom,Tv);


% Plotting - Total Variance
hold on;
h= errorbar(uts,vXv,seV,'k-v','LineWidth',3,'Markersize',12,'MarkerFaceColor','k'); % empirical data
plot(ts,V,'Color',[.4 .4 .4],'LineWidth',5); % DCDM variance
plot(ts,fvar(pvarcom,ts),'b','LineWidth',5); % quadratic variance for Gaussian models
plot([0,4],[6.1 6.1],'k-.','LineWidth',5,'Color',[.8 .8 .8]) % constant variance for Gaussian models
set(gca,'Fontsize',30);
axis([0 3.8 0 16])
legend('Empirical','DCDM','f(t)','Constant')

drawnow;
disp('Plotting of Variance Finished')


% Explicit variance function (see Supplemental Material)
function v = VarFun(r1,r2,x0,xeq,sd,Ts)
    A = r2/(r2-r1);
    a1 = 1/(r2-2*r1);
    a2 = 2/r1;
    a3 = -1/r2;

    E1 = a1*exp(-2*r1*Ts);
    E2 = a2*exp(-(r1+r2)*Ts);
    E3 = a3*exp(-2*r2*Ts);
    E4 = -(a1 + a2 + a3)*exp(-r2*Ts);

    v = sd^2 + r2*(x0-xeq)^2*(1-A)^2*(E1 + E2 + E3 + E4);
end