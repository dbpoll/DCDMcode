%% Generates Optimal Parameters from Table 1 
%
% Required Files:
% - face_dataset_paper.xlsx
% 
% Required Toolboxes:
% - Optimization Toolbox
% - Statistics and Machine Learning Toolbox
% - Financial Toolbox
% - Curve Fitting Toolbox
% 
% Code tested on MATLAB R2024b
       
%% Initialize Script
clear; clc; % clear command and figure windows

exporttable = false; % set whether table is exported to xlsx file

% Options for optimization and statistics toolboxes
opts1=  optimset('TolFun',1e-10,'TolX',1e-10,'display','off'); % least-squares display options
opts2 = statset('MaxIter',1e5,'MaxFunEvals',1e5,'Display','off'); % MLE options
opts2b = statset('MaxIter',1e5,'MaxFunEvals',1e5,'Display','off', ...
        'TolX',1e-12,'TolFun',1e-12,'TolTypeX','abs','TolTypeFun','abs'); % MLE options
opts3 = odeset('RelTol',1e-10,'AbsTol',1e-12); % ODE options


%% Import and Initialize Data
dtab = readtable('face_dataset_paper.csv'); % read data table
tshft = 18; % time offset from initial data collection and fumigation
t0 = min(dtab.DateSampled)-tshft; % set t = 0 as cessation of fumigation;
t = daysact(t0,dtab.DateSampled)/365; % set time (in years) for samples

% Indexing
Cind1 = strcmp(dtab.CO2Treatment,'Ambient'); % carbon treatment 
Cind2 = strcmp(dtab.CO2Treatment,'Elevated');

% Fixed Parameters 
T = 4.1; % end time (in years)
xeq = mean(dtab.d13C(Cind1));  % d13C equilibrium, estimated from ambient samples
n = 1; % number of digits for rounding in log-likelihood of table
nrnd = 10; % digits for rounding from data import

uTs = {3}; % unique time points of data collection
Ts = {3}; % initialize cell arrays for each order category, time
Xs = {3}; % initialize cell arrays for each order category, d13C
for i = 1:3
    if i == 1, Oind = (dtab.Order == 1) | (dtab.Order == 2); end % Orders 1 and 2
    if i == 2, Oind = (dtab.Order == 3) | (dtab.Order == 4); end % Orders 2 and 3
    if i == 3, Oind = (dtab.Order == 5); end % Order 5

    Ts{i} = t(Cind2 & Oind); % cell array for time data
    uTs{i} = unique(Ts{i}); % unique time points
    Xs{i} = round(dtab.d13C(Cind2 & Oind),nrnd); % cell array for d13C data
end

Tv = [Ts{1};Ts{2};Ts{3}]; % combined time data
Xv = [Xs{1};Xs{2};Xs{3}]; % combined d13C data

%% Variance Calculation from Data
vXv = zeros(1,9);
uts = unique(Tv);
for k = 1:length(uts)
    vXv(k) = var( Xv(Tv == uts(k) ) ); % compute variance each time point
end

% fit variance function to total empircal data
lbv = [-Inf -Inf -Inf];
ubv = [Inf Inf Inf];
f0 = [-2 10 0.5]; 
fcom = @(x) x(1).*uts.^2 + x(2).*uts + x(3) - vXv';

pvarcom = lsqnonlin(fcom,f0,lbv,ubv,opts1);

% Generate variance at time points
fvar = @(x,t) x(1).*t.^2 + x(2).*t + x(3);
varGauss = fvar(pvarcom,Tv);

disp('Variance for Data (Order Separated) Finished')


%% Simulation Comparisons 
% Assume:
% - xeq is known, fit from control (ambient) data
% - x0 is unknown
% - rates are unknown
% - standard deviation is unknown and is either constant or fit with a
%   low-order polynomial that maximizes the log-likelihood

% Single Rate Gaussian Fit (independent - r2)
gsi = @(x) [xeq + (x(end)-xeq)*exp(-x(1)*Ts{1}) - Xs{1}; ... % Orders 1 and 2
          xeq + (x(end)-xeq)*exp(-x(2)*Ts{2}) - Xs{2}; ... % Orders 3 and 4
          xeq + (x(end)-xeq)*exp(-x(3)*Ts{3}) - Xs{3}]; % Order 5
g0si = [.4, 0.4, 0.4, -38.4];
lbsi = [0 0 0 -50];
ubsi = [10 10 10 -32];
pfitsi = lsqnonlin(gsi,g0si,lbsi,ubsi,opts1);

% Maximum-likelihood 
psi = mle(Xv,'pdf',@(x,r_1,r_2,r_3,x0,sd)SEMi(x,r_1,r_2,r_3,x0,sd,xeq,Ts{1},Ts{2},Ts{3}),'Start',[pfitsi 2],'Options',opts2);
mSEMi = [xeq + (psi(end-1)-xeq)*exp(-psi(1)*Ts{1}); ... 
          xeq + (psi(end-1)-xeq)*exp(-psi(2)*Ts{2}); ... 
          xeq + (psi(end-1)-xeq)*exp(-psi(3)*Ts{3})]; 


SEMiLL = round( sum( log(normpdf(Xv,mSEMi,psi(end))) ),n );
SEMiAIC = 2*5- 2*SEMiLL; % total AIC

SEMiLLv = round(sum( log(normpdf(Xv,mSEMi,sqrt(varGauss))) ),n );
SEMiAICv = 2*(5-1+3) - 2*SEMiLLv; % total AIC

disp('Single Exponential Model (Independent Rates) Finished')


%% Single Rate Gaussian Fit (Dependent r1)
% Initialized with least-squares estimate
gsd = @(x) xeq + (x(end)-xeq)*exp(-x(1)*Tv) - Xv;
g0sd = [.4, -38.4];
lbsd = [0 -50];
ubsd = [10 -32];
pfitsd = lsqnonlin(gsd,g0sd,lbsd,ubsd,opts1); 

% Maximum-likelihood
psd = mle(Xv,'pdf',@(x,r,x0,sd)SEMd(x,r,x0,sd,xeq,Tv),...
    'Start',[pfitsd 2],'LowerBound',[lbsd 0],'UpperBound',[ubsd 100],'Options',opts2); 
mSEMd = xeq + (psd(2)-xeq)*( exp(-psd(1)*Tv) ); % mean of single-pool Gaussian

% Log-Likelihood and AIC
SEMdLL = round( sum( log(normpdf(Xv,mSEMd,psd(end))) ),n ); % log-likelihood of single-pool Gaussian
SEMdAIC = 2*3- 2*SEMdLL;  % total AIC

SEMdLLv = round( sum( log(normpdf(Xv,mSEMd,sqrt(varGauss))) ),n ); % log-likelihood of single-pool Gaussian
SEMdAICv = 2*(3-1+3)- 2*SEMdLLv;  % total AIC

disp('Single Exponential Model (Dependent Rates) Finished')


%%  Normally Distributed Double Exponential Model (Independent r1)
% Fitting strategy similar to Matamala or Lynch
% r1_1 - r1 for orders 1,2
% r1_2 - r1 for orders 3,4
% r1_3 - r1 for order 5
% r2_1 - r2 for orders 1,2
% r2_2 - r2 for orders 3,4
% r2_3 - r2 for order 5
% 
% a1_1 - coefficient for r1 for orders 1,2
% a1_2 - coefficient for r1 for orders 3,4
% a1_3 - coefficient for r1 for orders 5
% a2_1 - coefficient for r2 for orders 1,2
% a2_2 - coefficient for r2 for orders 3,4
% a2_3 - coefficient for r2 for orders 5

% x0 - initial mean delta C13, shared between datasets
% sd - Standard deviation of Gaussian distribution

gi = @(x) [xeq + (x(end)-xeq)*( x(7)*exp(-x(1)*Ts{1}) + x(10)*exp(-x(4)*Ts{1}) ) - Xs{1}; ... % Orders 1 and 2
          xeq + (x(end)-xeq)*( x(8)*exp(-x(2)*Ts{2}) + x(11)*exp(-x(5)*Ts{2}) ) - Xs{2};  % Orders 3 and 4
          xeq + (x(end)-xeq)*( x(9)*exp(-x(3)*Ts{3}) + x(12)*exp(-x(6)*Ts{3}) ) - Xs{3}]; % Order 5
g0i = [0.7, 0.7, 0.7, 3.6, 3.6, 3.6, 0.8, 0.8, 0.8, 0.2, 0.2, 0.2, -39];
lbi = [0 0 0 0 0 0 0 0 0 0 0 0 -50];
ubi = [10 10 10 100 10 10 10 10 10 10 10 10 -32];
pfitti = lsqnonlin(gi,g0i,lbi,ubi,opts1);

% Maximum-likelihood
pdi = mle(Xv,'pdf',@(x,r1_1,r1_2,r1_3,r2_1,r2_2,r2_3,a1_1,a1_2,a1_3,a2_1,a2_2,a2_3,x0,sd)...
        DEMi(x,r1_1,r1_2,r1_3,r2_1,r2_2,r2_3,a1_1,a1_2,a1_3,a2_1,a2_2,a2_3,x0,sd,xeq,Ts{1},Ts{2},Ts{3}),...
        'Start',[pfitti 2],'LowerBound',[lbi 0],'UpperBound',[ubi 100],'Options',opts2);

mDEMi = [xeq + (pdi(end-1)-xeq)*( pdi(7)*exp(-pdi(1)*Ts{1}) + pdi(10)*exp(-pdi(4)*Ts{1}) ); ... % Orders 1 and 2
    xeq + (pdi(end-1)-xeq)*( pdi(8)*exp(-pdi(2)*Ts{2}) + pdi(11)*exp(-pdi(5)*Ts{2}) ); ... % Orders 3 and 4
    xeq + (pdi(end-1)-xeq)*( pdi(9)*exp(-pdi(3)*Ts{3}) + pdi(12)*exp(-pdi(6)*Ts{3}) )]; % Order 5

% Log-likelihood and AIC
DEMiLL = round( sum( log(normpdf(Xv,mDEMi,pdi(end))) ), n);
DEMiAIC = 2*14- 2*DEMiLL; % total AIC

DEMiLLv = round( sum( log(normpdf(Xv,mDEMi,sqrt(varGauss))) ), n);
DEMiAICv = 2*(14-1+3)- 2*DEMiLLv; % total AIC

disp('Double Exponential Model (Independent Rates) Finished')


%% Normally Distributed Double Exponential Model (Dependent r1)
% Run Least-Squares for good initial conditions
% Separated by Order 1,2; 3,4; 5
% Parameters: [r1, r2_1, r2_2, r2_3, x0, sd]
% r1 - Carbon flux rate within tree - fixed across all orders
% r2_1 - Turnover rate for Order 1,2 roots
% r2_2 - Turnover rate for Order 3,4 roots
% r2_3 - Turnover rate for Order 5 roots
% x0 - initial mean delta C13
% xeq - equilibrium mean delta C13
% sd - Standard deviation of Gaussian distribution

g = @(x) [xeq + (x(end)-xeq)*( ( x(2)/(x(2) - x(1)) )*exp(-x(1)*Ts{1}) + ...
          ( 1 - x(2)/(x(2)-x(1)) )*exp(-x(2)*Ts{1}) ) - Xs{1}; ... % Orders 1 and 2
          xeq + (x(end)-xeq)*( ( x(3)/(x(3) - x(1)) )*exp(-x(1)*Ts{2}) + ...
          ( 1 - x(3)/(x(3)-x(1)) )*exp(-x(3)*Ts{2}) ) - Xs{2};  % Orders 3 and 4
          xeq + (x(end)-xeq)*( ( x(4)/(x(4) - x(1)) )*exp(-x(1)*Ts{3}) + ...
          ( 1 - x(4)/(x(4)-x(1)) )*exp(-x(4)*Ts{3}) ) - Xs{3}]; % Order 5
g0 = [0.4, 0.5, 0.5, 7, -38.4];
lb = [0 0 0 0 -50];
ub = [10 100 10 10 -32];
pfitt = lsqnonlin(g,g0,lb,ub,opts1);

% Maximum-likelihood 
pdd = mle(Xv,'pdf',@(x,r1,r2_1,r2_2,r2_3,x0,sd)DEMd(x,r1,r2_1,r2_2,r2_3,x0,sd,xeq,Ts{1},Ts{2},Ts{3}),...
        'Start',[pfitt 2],'LowerBound',[lb 0],'UpperBound',[ub 100],'Options',opts2);
a1 = pdd(2)/(pdd(2) - pdd(1));
a2 = pdd(3)/(pdd(3) - pdd(1));
a3 = pdd(4)/(pdd(4) - pdd(1));
mDEMd = [xeq + (pdd(end-1)-xeq)*( a1*exp(-pdd(1)*Ts{1}) + (1-a1)*exp(-pdd(2)*Ts{1}) ); ... % Orders 1 and 2
    xeq + (pdd(end-1)-xeq)*( a2*exp(-pdd(1)*Ts{2}) + (1-a2)*exp(-pdd(3)*Ts{2}) ); ... % Orders 3 and 4
    xeq + (pdd(end-1)-xeq)*( a3*exp(-pdd(1)*Ts{3}) + (1-a3)*exp(-pdd(4)*Ts{3}) )]; % Order 5

% Log-likelihood and AIC
DEMdLL = round( sum( log(normpdf(Xv,mDEMd,pdd(end))) ), n); 
DEMdAIC = 2*6- 2*DEMdLL; % total AIC

DEMdLLv = round( sum( log(normpdf(Xv,mDEMd,sqrt(varGauss))) ), n); 
DEMdAICv = 2*(6-1+3)- 2*DEMdLLv; % total AIC

disp('Double Exponential Model (Dependent Rates) Finished')

%% Demographic Carbon Distribution Model
% Use Least-Squares from Two Pool Gaussian as initial conditions
% Separated by Order 1,2; 3,4; 5
% Parameters: pdcdm = [r1, r2_1, r2_2, r2_3, x0, sd]
% r1 - Carbon flux rate within tree - fixed across all orders
% r2_1 - Turnover rate for Order 1,2 roots
% r2_2 - Turnover rate for Order 3,4 roots
% r2_3 - Turnover rate for Order 5 roots
% x0 - initial mean delta C13
% xeq - equilibrium mean delta C13

% Maximum Likelihood
p_temp = pdd; p_temp(2) = 1;
pdcdm = mle(Xv,'pdf',@(x,r1,r2_1,r2_2,r2_3,x0,sd)DCDM(x,r1,r2_1,r2_2,r2_3,x0,sd,xeq,Ts{1},Ts{2},Ts{3},opts3),...
    'Start',p_temp,'LowerBound',[lb 0],'UpperBound',[ub 100],'Options',opts2);

% Log-likelihood and AIC
DCDMLL = 0;
for i = 1:3
    DCDMLL = DCDMLL + DCDMfunLL(Xs{i},pdcdm(1),pdcdm(i+1),pdcdm(end),pdcdm(end-1),xeq,Ts{i},opts3);
end
DCDMLL = round( DCDMLL, n);
DCDMAIC = 2*6- 2*DCDMLL; % total AIC
disp('DCDM Optimization Finished')



%% LOESS Regression
plo = mle(Xv,'pdf',@(x,b_1,b_2,b_3,sd)LoMod(x,b_1,b_2,b_3,sd,Ts{1},Ts{2},Ts{3}),...
         'Start',[0.2 0.2 0.2 1.5],...
        'LowerBound',[0 0 0 0.1],...
        'UpperBound',[1 1 1 Inf],'Options',opts2b);

LoLL = round( sum( log(LoMod(Xv,plo(1),plo(2),plo(3),plo(end),Ts{1},Ts{2},Ts{3})) ),n);
LoAIC = 2*(3*6) - 2*LoLL;

LoLLv = round( sum( log(LoMod(Xv,plo(1),plo(2),plo(3),sqrt(varGauss),Ts{1},Ts{2},Ts{3})) ),n );
LoAICv = 2*(3*6-1+4) - 2*LoLLv;

disp('LOESS Model Finished')


%% Construct Figure Table
Model = ["SEM_dependent";"";"SEM_independent";"";"DEM_dependent";"";"DEM_independent";""; ...
         "LOESS";"";"DCDM"];
Order12_r1 = [psd(1); NaN; psi(1); NaN; pdd(1); NaN; pdi(1); NaN; NaN; NaN; pdcdm(1)];
Order12_r2 = [NaN; NaN; NaN; NaN; pdd(2); NaN; pdi(4); NaN; NaN; NaN; pdcdm(2)];
Order34_r1 = [psd(1); NaN; psi(2); NaN; pdd(1); NaN; pdi(2); NaN; NaN; NaN; pdcdm(1)];
Order34_r2 = [NaN; NaN; NaN; NaN; pdd(3); NaN; pdi(5); NaN; NaN; NaN; pdcdm(3)];
Order5_r1 = [psd(1); NaN; psi(3); NaN; pdd(1); NaN; pdi(3); NaN; NaN; NaN; pdcdm(1)];
Order5_r2 = [NaN; NaN; NaN; NaN; pdd(4); NaN; pdi(6); NaN; NaN; NaN; pdcdm(4)];
x_0 = [psd(2); NaN; psi(4); NaN; pdd(5); NaN; pdi(13); NaN; NaN; NaN; pdcdm(5)];
sigma2 = [psd(3).^2; NaN; psi(5).^2; NaN; pdd(6).^2; NaN; pdi(14).^2; NaN; NaN; NaN; pdcdm(6).^2];
LL = [SEMdLL; SEMdLLv; SEMiLL; SEMiLLv; DEMdLL; DEMdLLv; DEMiLL; DEMiLLv; LoLL; LoLLv; DCDMLL];
AIC = [SEMdAIC; SEMdAICv; SEMiAIC; SEMiAICv; DEMdAIC; DEMdAICv; DEMiAIC; DEMiAICv; LoAIC; LoAICv; DCDMAIC];
FP = [3; 5; 5; 7; 6; 8; 14; 16; 18; 20; 6];

table = table(Model,Order12_r1,Order12_r2,Order34_r1,Order34_r2,Order5_r1,Order5_r2,...
            x_0,sigma2,LL,AIC,FP);

disp(table)
if exporttable
    writetable(table,'table1.xlsx');
    disp('Table written to table1.xlsx');
end



%% Auxillary Functions
% Functions are:
% - Single-Pool Gaussian Model (independent) 
% - Single-Pool Gaussian Model (dependent) 
% - Two-Pool Gaussian Model (independent) 
% - Two-Pool Gaussian Model (dependent) 
% - Demographic Carbon Distribution Model (dependent) 
% - LOESS Model (independent)

%% Single Pool Independent Functions
% Compute mean, then run data through normal PDF
function SEMi = SEMi(x,r_1,r_2,r_3,x0,sd,xeq,Ts_1,Ts_2,Ts_3)
    mSEM =  [xeq + (x0-xeq)*exp(-r_1*Ts_1); ... % Orders 1 and 2
          xeq + (x0-xeq)*exp(-r_2*Ts_2); ... % Orders 3 and 4
          xeq + (x0-xeq)*exp(-r_3*Ts_3)]; % Order 5
   SEMi = normpdf(x,mSEM,sd);
end


%% Single Pool Dependent Functions 
% Compute mean, then run data through normal PDF
function SEMd = SEMd(x,r,x0,sd,xeq,Tv)
    mSEM = xeq + (x0-xeq)*( exp(-r*Tv) ); % mean of two-pool
    SEMd = normpdf(x,mSEM,sd);
end


%% Two Pool Gaussian (Independent) Functions
% Two Pool Model Function
function DEMi = DEMi(x,r1_1,r1_2,r1_3,r2_1,r2_2,r2_3,a1_1,a1_2,a1_3,a2_1,a2_2,a2_3,x0,sd,xeq,Ts_1,Ts_2,Ts_3)

    mDEM = [xeq + (x0-xeq)*( a1_1*exp(-r1_1*Ts_1) + a2_1*exp(-r2_1*Ts_1) ); ... % Orders 1 and 2
          xeq + (x0-xeq)*( a1_2*exp(-r1_2*Ts_2) + a2_2*exp(-r2_2*Ts_2) ); ... % Orders 3 and 4
          xeq + (x0-xeq)*( a1_3*exp(-r1_3*Ts_3) + a2_3*exp(-r2_3*Ts_3) )]; % Order 5

    DEMi=normpdf(x,mDEM,sd);
end

%% Two Pool Gaussian (Dependent) Functions
% Two Pool Model Function
function DEMd = DEMd(x,r1,r2_1,r2_2,r2_3,x0,sd,xeq,Ts_1,Ts_2,Ts_3)
    a1 = r2_1/(r2_1 - r1);
    a2 = r2_2/(r2_2 - r1);
    a3 = r2_3/(r2_3 - r1);
    mDEM = [xeq + (x0-xeq)*( a1*exp(-r1*Ts_1) + (1-a1)*exp(-r2_1*Ts_1) ); ... % Orders 1 and 2
          xeq + (x0-xeq)*( a2*exp(-r1*Ts_2) + (1-a2)*exp(-r2_2*Ts_2) ); ... % Orders 3 and 4
          xeq + (x0-xeq)*( a3*exp(-r1*Ts_3) + (1-a3)*exp(-r2_3*Ts_3) )]; % Order 5

    DEMd=normpdf(x,mDEM,sd);
end


%% Two Pool Markov (Our Model) Function
function DCDM = DCDM(x,r1,r2_1,r2_2,r2_3,x0,sd,xeq,Ts_1,Ts_2,Ts_3,opts)
    x = reshape(x,[],1); % forces column vector
    T = 1500;

    ind1 = length(Ts_1);
    ind2 = length(Ts_2);
    ind3 = length(Ts_3);
    DCDM = [];
    for i = 1:3 % get r2, Ts{i}, and Xs{i}
        if i == 1, r2 = r2_1; to = Ts_1; xo = x(1:ind1); end
        if i == 2, r2 = r2_2; to = Ts_2; xo = x((ind1+1):(ind1+ind2)); end
        if i == 3, r2 = r2_3; to = Ts_3; xo = x((ind1+ind2+1):(ind1+ind2+ind3)); end
        v = [r1 r2 sd];

        P0 = normpdf(xo,x0,sd); % define initial distribution
        [s,P] = ode45(@(t,P) Pfun(xo,t,P,v,x0,xeq),[0 T],P0,opts); % run ODE solver
        Pint = diag(interp1(s,P,to,'spline')); % interpolate to time-points
        DCDM = [DCDM;Pint];
    end

end

% Diff Eq. for Our Model
function dPdt = Pfun(x,t,P,v,x0,xeq)
    x = reshape(x,[],1); % forces column vector    
    r1 = v(1);
    r2 = v(2);
    sd = v(3);

    B = xeq + (x0-xeq)*exp(-r1*t);
    dPdt = -r2*P + r2*normpdf(x,B,sd);

end

% Log-likelihood Function Calculation
function DCDMfunLL = DCDMfunLL(x,r1,r2,sd,x0,xeq,ts,opts)
    x = reshape(x,[],1); % forces column vector
    T = 1500; 
    v = [r1 r2 sd];
    P0 = normpdf(x,x0,sd);
    [s,P] = ode45(@(t,P) Pfun(x,t,P,v,x0,xeq),[0 T],P0,opts);
    Pint = interp1(s,P,ts,'spline');
    DCDMfunLL = sum( log(diag(Pint)) );
end



%% Loess Model Functions
% Run data through normal PDF
function Lo = LoMod(x,b_1,b_2,b_3,sd,Ts_1,Ts_2,Ts_3)
    ind1 = length(Ts_1);
    ind2 = length(Ts_2);
    ind3 = length(Ts_3);
    z = zeros(length(x),1);
    for i = 1:3
        if i == 1, z(1:ind1) = smooth(Ts_1,x(1:ind1),b_1,'loess'); end
        if i == 2, z((ind1+1):(ind1+ind2)) = smooth(Ts_2,x((ind1+1):(ind1+ind2)),b_2,'loess'); end
        if i == 3, z((ind1+ind2+1):(ind1+ind2+ind3)) = smooth(Ts_3,x((ind1+ind2+1):(ind1+ind2+ind3)),b_3,'loess');  end
    end
    Lo = normpdf(x,z,sd);
end