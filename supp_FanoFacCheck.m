% Run Fano factor check to ID how far each cortical area departs from the Poisson assumptions under 
% our FI calculation and make sure there aren't substantial difference from
% each other

clear all; close all;

% Load in processed data
load('./analysisFiles/physio/fittingResults_processed.mat');

%% Concatenate data 
% (mean & variance for each disparity presented to each cell - 
% (i.e. one cell's data appears more than once)

% initialize data structures to store means/variances/repeats
V1dat = []; V2dat = []; MTdat = [];
MTids = []; V1ids = []; V2ids = [];

% grab means (col2), variances (col3), and repeats (col4) for each cell in each area
for ii = 1:numel(V1.experiments)
    V1dat = [V1dat; V1.experiments{ii}.dat(:,2:4)];
    V1ids = [V1ids; ii*ones(size(V1.experiments{ii}.dat(:,2:4),1),1)];
end

for ii = 1:numel(V2.experiments)
    V2dat = [V2dat; V2.experiments{ii}.dat(:,2:4)];
    V2ids = [V2ids; ii*ones(size(V2.experiments{ii}.dat(:,2:4),1),1)];
end

for ii = 1:numel(MT.experiments)
    MTdat = [MTdat; MT.experiments{ii}.dat(:,2:4)];
    MTids = [MTids; ii*ones(size(MT.experiments{ii}.dat(:,2:4),1),1)];
end

%% Cull data based on number of repeats and maximum variance
minReps = 3;

V1dat = V1dat(V1dat(:,3) >= minReps,:);
V2dat = V2dat(V2dat(:,3) >= minReps,:);
MTdat = MTdat(MTdat(:,3) >= minReps,:);

V1ids = V1ids(V1dat(:,3) >= minReps,:);
V2ids = V2ids(V2dat(:,3) >= minReps,:);
MTids = MTids(MTdat(:,3) >= minReps,:);

% replace zeros with eps to allow fitting of power function
V1dat(V1dat(:,2) == 0,2) = eps;
V2dat(V2dat(:,2) == 0,2) = eps;
MTdat(MTdat(:,2) == 0,2) = eps;

V1dat(V1dat(:,1) == 0,1) = eps;
V2dat(V2dat(:,1) == 0,1) = eps;
MTdat(MTdat(:,1) == 0,1) = eps;

% Find number of cells from each area
MTcells = numel(unique(MTids));
V1cells = numel(unique(V1ids));
V2cells = numel(unique(V2ids));


%% Find best fitting power law

% specify fitting routine
pLaw   = @(a,x,b) a*x.^b;
p0     = [1 1];
opts   = optimoptions('fmincon','Display','off');
A      = [0,-1];
b      = -1;

% concatenate the data
catDat = {V1dat,V2dat,MTdat};

% for each area
for ii = 1:3

    % grab the data from that area
    data = catDat{ii};

    % specify objective function
    objFxn = @(p) sum( ( log(pLaw(p(1),data(:,1),p(2))) - log(data(:,2)) ).^2 );
    
    % run fitting
    pFit(ii,:) = fmincon(objFxn,p0,A,b,[],[],[],[],[],opts);
end


%% Plot

% support for plotting
supp = linspace(1e-1,1e4,100);

xtix = [0.01 0.1 1 10 100 1000];
for ii = 1:numel(xtix)
    xLab{ii} = num2str(xtix(ii));
end

f1 = figure;
f1.Position = [2000 100 1700 450];

% V1 distribution
subplot(1,3,1);
hold on;

scatter(V1dat(:,1),V1dat(:,2),70,[0.6 0.6 0.6]);
plot(supp,pLaw(pFit(1,1),supp,pFit(1,2)),'r','LineWidth',2);
plot([1e-1 10000],[1e-1 10000],'--k','LineWidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLab,'xlim',[1e-1 10000],'ylim',[1e-1 10000],'fontsize',15,'yscale','log','xscale','log');
xlabel('mean spike count');
ylabel('spike count variance');
title(['V1 (n=',num2str(V1cells),')']);

% V2 distribution
subplot(1,3,2);
hold on;

scatter(V2dat(:,1),V2dat(:,2),70,[0.6 0.6 0.6]);
plot(supp,pLaw(pFit(2,1),supp,pFit(2,2)),'r','LineWidth',2);
plot([1e-1 10000],[1e-1 10000],'--k','LineWidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLab,'xlim',[1e-1 10000],'ylim',[1e-1 10000],'fontsize',15,'yscale','log','xscale','log');
xlabel('mean spike count');
ylabel('spike count variance');
title(['V2 (n=',num2str(V2cells),')']);

% MT distribution
subplot(1,3,3);
hold on;

scatter(MTdat(:,1),MTdat(:,2),70,[0.6 0.6 0.6]);
plot(supp,pLaw(pFit(3,1),supp,pFit(3,2)),'r','LineWidth',2);
plot([1e-1 10000],[1e-1 10000],'--k','LineWidth',2);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',xtix,'xticklabel',xLab,'xlim',[1e-1 10000],'ylim',[1e-1 10000],'fontsize',15,'yscale','log','xscale','log');
xlabel('mean spike count');
ylabel('spike count variance');
title(['MT (n=',num2str(MTcells),')']);

% Plot best-fitting power law
f2 = figure;
f2.Position = [2000 400 1200 450];

subplot(1,2,1);
hold on;
X = categorical({'V1','V2','MT'});
X = reordercats(X,{'V1','V2','MT'});
Y = pFit(:,1);
bar(X,Y,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
plot([0 4],[1 1],'--k','linewidth',2);
ylabel('Slope');
set(gca,'plotboxaspectratio',[1 1 1],'xlim',{'V1','MT'},'fontsize',15);

subplot(1,2,2);
hold on;
Y = pFit(:,2);
bar(X,Y,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
plot([0 4],[1 1],'--k','linewidth',2);
ylabel('Power law');
set(gca,'plotboxaspectratio',[1 1 1],'xlim',{'V1','MT'},'fontsize',15);

%% Save figures

saveas(f1,'./plots/FanoFactorCheck/FanoFactor_scatter.svg');
saveas(f2,'./plots/FanoFactorCheck/FanoFactor_fits.svg');

