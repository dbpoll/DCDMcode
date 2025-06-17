%% Figure 1 - Plot of Data and Variance%
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

clc; clf; % clear command and figure windows
clearvars -except dtab; % clear variables except table
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');  % removes table header warning

opts1=  optimset('TolFun',1e-10,'TolX',1e-10,'display','off'); % least-squares display options
opts2 = statset('MaxIter',1e5,'MaxFunEvals',1e5,'Display','off'); % MLE options
opts3 = odeset('RelTol',1e-10,'AbsTol',1e-12); % ODE options

%% Import and Initialize Data
dtab = readtable('face_dataset_paper.xlsx'); % read d13C data table 
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

Ts = {3};
Xs = {3};
for i = 1:3
    % Categorize data into Orders [1,2] [3,4] and [5]
    if i == 1, Oind = (dtab.Order == 1) | (dtab.Order == 2); end % Orders 1 and 2
    if i == 2, Oind = (dtab.Order == 3) | (dtab.Order == 4); end % Orders 2 and 3 
    if i == 3, Oind = (dtab.Order == 5); end % Order 5

    Ts{i} = t(Cind2 & Oind);
    Xs{i} = dtab.d13C(Cind2 & Oind);

end

%% Variance Calculation from Data
vXsT = zeros(3,9);
seT = zeros(3,9);
for i = 1:3
    % Variance Calculation - Separated by Order/Diam
    uts = unique(Ts{i}); % find unique time values
    vXs = zeros(1,length(uts));
    cnt_pts = zeros(1,length(uts));
    mu4 = zeros(1,length(uts));
    for k = 1:length(uts)
        cnt_pts(k) =  length(Xs{i}(Ts{i} == uts(k)));
        vXs(k) = var( Xs{i}(Ts{i} == uts(k) ) ); % variance at unique time values
        mu4(k) = mean( ( Xs{i}(Ts{i} == uts(k)) - mean( Xs{i}(Ts{i} == uts(k)) )  ).^4 ); % mu4 from Rao 73
    end
    vXsT(i,:) = vXs;
    seT(i,:) = sqrt(  (mu4 - (cnt_pts-3).*(vXs.^2) ./(cnt_pts-1) )./(cnt_pts)  );
end
disp('Variance for Data (Order Separated) Finished')

% Total Variance Calculuation
Tv = [Ts{1};Ts{2};Ts{3}]; % combine time arrays
Xv = [Xs{1};Xs{2};Xs{3}]; % combine d13C arrays

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

%% Variance Calculation from Model

% Run Model - Separated by Order [1,2], [3,4], and [5]
% Parameters: [r1, r2_1, r2_2, r2_3, x0, sd]
% r1 - Effective rate of isotopes in tree carb pool
% r2_1 - Turnover rate for order 1,2 roots
% r2_2 - Turnover rate for order 3,4 roots
% r2_3 - Turnover rate for order 5 roots
% x0 - initial mean delta C13
% xeq - equilibrium mean delta C13

% % Maximum Likelihood
p0 = [0.56 1 4.4 9 -38.4 2.5];
lb = [0 0 0 0 -50 0];
ub = [10 100 10 10 -32 100];
pDCDM = mle(Xv,'pdf',@(x,r1,r2_1,r2_2,r2_3,x0,sd)DCDM(x,r1,r2_1,r2_2,r2_3,x0,sd,xeq,Ts{1},Ts{2},Ts{3},opts3),...
    'Start',p0,'LowerBound',lb,'UpperBound',ub,'Options',opts2);

% Calculate Variance from model
v = zeros(3,length(ts));
for i = 1:3
    v(i,:) = VarFun(pDCDM(1),pDCDM(i+1),pDCDM(end-1),xeq,pDCDM(end),ts);
end

disp('Variance for DCDM (Order Separated) Finished')

% Run Model - Total (Combined)
% Parameters: [r1, r2, x0, sd]
% r1 - Effective rate of isotopes in tree carb pool
% r2 - Average turnover rate for all orders 1-5
% x0 - initial mean delta C13
% xeq - equilibrium mean delta C13

q0 = [pDCDM(1) mean(pDCDM(2:4)) pDCDM(5)  pDCDM(6)];
lbc = [0 0 -50 0];
ubc = [10 100 -32 100];
pDCDMcom = mle(Xv,'pdf',@(x,r1,r2,x0,sd)DCDMcom(x,r1,r2,x0,sd,xeq,Tv,opts3),...
    'Start',q0,'LowerBound',lbc,'UpperBound',ubc,'Options',opts2);


% Calculate Variance
V = VarFun(pDCDMcom(1),pDCDMcom(2),pDCDMcom(end-1),xeq,pDCDMcom(end),ts);

disp('Variance for DCDM (Total) Finished')

%% Plotting Figure Panels
% Define figure parameters
cmap = cool(3); % array of color values for data plot
cmap2 = [0 0 255; 76 0 153; 204 0 0]/255; % color for model lines

for i =1:3
    % Plotting - Data (Panel C)
    subplot(3,2,[3 5])
    hold on;
    plot(Ts{i},Xs{i},'o','Markersize',14,'Color',[0 0 0],'MarkerFaceColor',cmap(i,:))

    % Plotting = Variance (Panel D)
    subplot(3,2,2*i) % variance, separated by order
    hold on;
    errorbar(uts,vXs,seT(i,:),'-v','LineWidth',4,'Markersize',8,'Color',cmap(i,:),'MarkerFaceColor',cmap(i,:));
    plot(ts,v(i,:),'LineWidth',6,'Color',cmap2(i,:));

    axis([0 4 0 20])
    set(gca,'Fontsize',30);
end

% Plotting - Equilibrium (Panel C)
subplot(3,2,[3 5])
plot([0 T],[xeq xeq],'k--','LineWidth',4); % equilibrum
axis([0 4 -45 -25])
set(gca,'Fontsize',30);
legend('Orders 1,2','Orders 3,4','Order 5')

% Plotting - Total Variance (Panel B)
subplot(3,2,1) % total variance
hold on;
h= errorbar(uts,vXv,seV,'k-v','LineWidth',3,'Markersize',12,'MarkerFaceColor','k');
plot(ts,V,'Color',[.4 .4 .4],'LineWidth',5);
plot([0,4],[6.1 6.1],'k-.','LineWidth',5,'Color',[.8 .8 .8])
set(gca,'Fontsize',30);
axis([0 4 0 20])

drawnow;
disp('Plotting of Variance Finished')


%% Auxillary Functions
% Function to return DCDM Distribution - Order Separated 
function DCDM = DCDM(x,r1,r2_1,r2_2,r2_3,x0,sd,xeq,Ts_1,Ts_2,Ts_3,opts)
    x = reshape(x,[],1); % forces column vector
    T = 1500;

    ind1 = length(Ts_1);
    ind2 = length(Ts_2);
    ind3 = length(Ts_3);
    inds = [0,ind1,ind1+ind2,ind1+ind2+ind3];
    DCDM = zeros(ind1+ind2+ind3,1);
    for i = 1:3 % get r2, Ts{i}, and Xs{i}
        if i == 1, r2 = r2_1; to = Ts_1; end
        if i == 2, r2 = r2_2; to = Ts_2; end
        if i == 3, r2 = r2_3; to = Ts_3; end
        xo = x((inds(i)+1):inds(i+1));
        v = [r1 r2 sd];

        P0 = normpdf(xo,x0,sd); % define initial distribution
        [s,P] = ode45(@(t,P) Pfun(xo,t,P,v,x0,xeq),[0 T],P0,opts); % run ODE solver
        Pint = diag(interp1(s,P,to,'spline')); % interpolate to time-points
        DCDM((inds(i)+1):inds(i+1)) = Pint;
    end
end

% Function to return DCDM Distribution - Total
function DCDM = DCDMcom(x,r1,r2,x0,sd,xeq,Tv,opts)
    x = reshape(x,[],1); % forces column vector
    T = 1500; % end time
    v = [r1 r2 sd];

    P0 = normpdf(x,x0,sd); % define initial distribution
    [s,P] = ode45(@(t,P) Pfun(x,t,P,v,x0,xeq),[0 T],P0,opts); % run ODE solver
    DCDM = diag(interp1(s,P,Tv,'spline')); % interpolate to time-points

end

% Solve DCDM diff. equation
function dPdt = Pfun(x,t,P,v,x0,xeq)
    x = reshape(x,[],1); % forces column vector    
    r1 = v(1);
    r2 = v(2);
    sd = v(3);
    
    C = xeq + (x0-xeq)*exp(-r1*t);
    dPdt = -r2*P + r2*normpdf(x,C,sd);

end

% Explicit variance function (see App. B)
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
