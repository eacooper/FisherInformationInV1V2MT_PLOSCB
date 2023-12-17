clear all; close all;

%% Figure 1

%% 1A (bottom)
fig1A = figure;
fig1A.Position = [100 100 650 600];
hold on;

x     = linspace(-2,2,100);
tcMus = linspace(-1.65,1.65,5);
sig   = 0.4;

% Plot 5 example Gaussian tuning functions
for ii = 1:5

    tcTint = [1 1 1]*0.75*(ii-1)/5;
    tc = normpdf(x,tcMus(ii),sig);

    plot(x,tc,'color',tcTint,'linewidth',8);

end

set(gca,'xlim',[-2 2],'xtick',-2:2,'ylim',[0 1.1],'ytick',[],'fontsize',30,'plotboxaspectratio',[1 1 1]);
xlabel('Binocular disparity (\circ)');
ylabel('Firing rate');


%% 1C

res = 52; lb  = -2; ub  = 2;edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

s = 0.3;
p = 1;

p1 = exp( -(abs(cntr_disp)/s).^p );
p1 = p1 / sum(p1);

p2 = p1.^2;
p2 = p2 / sum(p2);

ppt5 = p1.^0.5;
ppt5 = ppt5 / sum(ppt5);

% Plot probabilities to 0.5 and 2 powers
fig1C = figure;
fig1C.Position = [800 100 650 600];
hold on;

plot(cntr_disp, p1,'-','color',[0 0 0],'linewidth',2);
plot(cntr_disp, p2,'--','color',[0 0 0],'linewidth',2);
plot(cntr_disp, ppt5,':','color',[0 0 0],'linewidth',2);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:2,'fontsize',30)
axis square; box on;
ylim([0 .27])
legend('p=1','infomax','discrimax');
xlabel('Binocular disparity (\circ)');
ylabel('Fisher information (FI)');

%% 1D

dispTix = [-2 0 2];

s = RandStream('mt19937ar','Seed',94720);
RandStream.setGlobalStream(s);

% Define simulated data points
numDat = 80;

xVals = linspace(-1.8,1.8,numDat);

% Set total reconstruction error x penalty
totalPenalty = 42;

% Generate some random data for each error function
%--------------------------------------
% Uniformly sample some reconstruction error
error = 4*rand([1 numDat]) - 2;

l0penalty = sum(error~=0);
l2penalty = sum(error.^2);

% l0: select some random inds to make perfect to match designated optimum
l0error              = error;
zeroErrInds          = datasample(1:numDat,l0penalty-totalPenalty,'replace',false);
l0error(zeroErrInds) = 0;

% l2: scale the input error to reach designated optimum
l2error              = error*sqrt(totalPenalty/(sum(error.^2)));

% Generate reconstructed datasets under each optimization routine (norm)
l0dat = xVals + l0error;
l2dat = xVals + l2error;

% L0 norm
fig1Di = figure;
fig1Di.Position = [300 300 650 600];
hold on;

plot(x,ones(1,100),'k','linewidth',4);
plot([0 0],[0 1],'k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[0 1],'ytick',[0 1],'fontsize',30);
ylabel('Error Penalty'); 
xlabel('Reconstruction error');

% L0 reconstruction
fig1Dii = figure;
fig1Dii.Position = [980 300 650 600];
hold on;

scatter(xVals,l0dat,200,'k','filled');
plot([-4 4],[-4 4],'--k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-4 4],'xtick',dispTix,'ylim',[-4 4],'ytick',dispTix,'fontsize',30);
ylabel('Reconstruction'); 
xlabel('Binocular disparity (\circ)');

% L0 reconstruction error
fig1Diii = figure;
fig1Diii.Position = [980 300 650 600];
hold on;

scatter(xVals,l0dat-xVals,200,'k','filled');
plot([-4 4],[0 0],'--k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[-4 4],'ytick',dispTix,'fontsize',30);
ylabel('Reconstruction Error'); 
xlabel('Binocular disparity (\circ)');

% L2 norm
fig1Div = figure;
fig1Div.Position = [300 800 650 600];
hold on;

plot(x,(x.^2)/max(x.^2),'k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[0 1],'ytick',[0 1],'fontsize',30);
ylabel('Error Penalty'); 
xlabel('Reconstruction error');

% L2 reconstruction
fig1Dv = figure;
fig1Dv.Position = [980 800 650 600];
hold on;

scatter(xVals,l2dat,200,'k','filled');
plot([-4 4],[-4 4],'--k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-4 4],'xtick',dispTix,'ylim',[-4 4],'ytick',dispTix,'fontsize',30);
ylabel('Reconstruction'); 
xlabel('Binocular disparity (\circ)');

% L2 reconstruction error
fig1Dvi = figure;
fig1Dvi.Position = [980 300 650 600];
hold on;

scatter(xVals,l2dat-xVals,200,'k','filled');
plot([-4 4],[0 0],'--k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[-4 4],'ytick',dispTix,'fontsize',30);
ylabel('Reconstruction Error'); 
xlabel('Binocular disparity (\circ)');

% Save figs

saveas(fig1A,['./plots/Figure1Panels/fig1A.svg']);
saveas(fig1C,['./plots/Figure1Panels/fig1C.svg']);

saveas(fig1Di,['./plots/Figure1Panels/fig1Di.svg']);
saveas(fig1Dii,['./plots/Figure1Panels/fig1Dii.svg']);
saveas(fig1Diii,['./plots/Figure1Panels/fig1Diii.svg']);
saveas(fig1Div,['./plots/Figure1Panels/fig1Div.svg']);
saveas(fig1Dv,['./plots/Figure1Panels/fig1Dv.svg']);
saveas(fig1Dvi,['./plots/Figure1Panels/fig1Dvi.svg']);


%% Figure 5A

clear all; close all;

% Generate an example Gabor tuning curve + component functions

clear all
close all

% Define Gabor tuning function and component functions
TC = @(disp,r0,A,prPosDisp,sig,prFreq,prPhase) ...
    r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 ).*...
    cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

Gaus = @(disp,r0,A,prPosDisp,sig) ...
    r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 );

cosF = @(disp,r0,A,prPosDisp,prFreq,prPhase) ...
    r0 + A*cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

% Define example tuning curve parameters
d   = linspace(-2,2,1000);
r0  = 30;
A   = 40;
pd  = 0;
sig = 1;
pf  = 0.2;
pPh = 0.8;

% Generate tuning curves
gabFxn = TC(d,r0,A,pd,sig,pf,pPh);
gauFxn = Gaus(d,0,1,pd,sig);
cosFxn = cosF(d,0,1,pd,pf,pPh);


% Plot
ymax = 80;

fig1 = figure;
fig1.Position = [100 100 650 600];

hold on;

plot(d,gabFxn,'k','linewidth',4);
plot([-2 2],r0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 ymax],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig2 = figure;
fig2.Position = [100 100 650 600];

hold on;

plot(d,cosFxn,'--k','linewidth',4);
plot(-pPh/(2*pi*pf)*[1 1],[-1 1],'--k','linewidth',2);
plot(0*[1 1],[-1 1],'--k','linewidth',2);
plot([-2 2],0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[-1 1],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig3 = figure;
fig3.Position = [100 100 650 600];

hold on;

plot(d,gauFxn,'--k','linewidth',4);
plot(pd*[1 1],[0 1],'--k','linewidth',2);
plot(sig*[1 1],[0 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 1],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');


% Save plots
saveas(fig1,['./plots/Figure5Apanel/fullTC.svg']);
saveas(fig2,['./plots/Figure5Apanel/cosTC.svg']);
saveas(fig3,['./plots/Figure5Apanel/gauTC.svg']);


%% Figure 6A

% Generate a pair of example Gabor tuning curves (original and resampled)

clear all; close all;

% Define Gabor tuning function
TC = @(disp,r0,A,prPosDisp,sig,prFreq,prPhase) ...
    r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 ).*...
    cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

% Original V1 Cell
d   = linspace(-2,2,1000);
r0  = 30;
A   = 30;
pd  = 0;
sig = 0.75;
pf = 0.4;
pPh = 0.8;

% Resampled V1 Cell (only envelope std from MT distribution)
sig2 = 1.5;

% Define tuning functions
gabFxn1 = TC(d,r0,A,pd,sig,pf,pPh);
gabFxn2 = TC(d,r0,A,pd,sig2,pf,pPh);

% Plot
ymax = 80;

fig1 = figure;
fig1.Position = [100 500 650 600];

hold on;

plot(d,gabFxn1,'k','linewidth',4);

set(gca,'xlim',[-2 2],'ylim',[0 ymax],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig2 = figure;
fig2.Position = [800 500 650 600];

hold on;

plot(d,gabFxn2,'k','linewidth',4);

set(gca,'xlim',[-2 2],'ylim',[0 ymax],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');


% Save plots
saveas(fig1,['./plots/Figure6Apanel/originalTC.svg']);
saveas(fig2,['./plots/Figure6Apanel/resampledTC.svg']);

