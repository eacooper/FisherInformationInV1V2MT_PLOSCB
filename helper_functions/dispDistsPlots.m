function [] = dispDistsPlots()
% Plot disparity probability distributions used for main figs in paper and controls

% Load in disparity histograms based on each KSD
load('./analysisFiles/disparityStats/dispHistStruct.mat');

% Define stats subsets
imgSet      = {'sando','walking'};
region      = {'Circ','V1','V2','MT'};

% Define some more pleasant colors
blueCI  = ColorIt('b');
greenCI = ColorIt('g');
redCI   = ColorIt('r');

f = cell(4,1);
fCirc = cell(2,1);

% Loop over behavioral image sets
for ii = 1:2

    % Plot stats using full MT dataset
    f{ii} = figure;
    f{ii}.Position = [100 700 650 600];
    hold on;

    plot([0 0],[0 5],'--k','linewidth',2);

    p(1) = shadeplot(realDat.(imgSet{ii}).(region{2}).('all'),[dispDat.(imgSet{ii}).(region{2}).('LCI'); ...
        dispDat.(imgSet{ii}).(region{2}).('UCI')],cntr_disp,blueCI,0.5,4);
    p(2) = shadeplot(realDat.(imgSet{ii}).(region{3}).('all'),[dispDat.(imgSet{ii}).(region{3}).('LCI'); ...
        dispDat.(imgSet{ii}).(region{3}).('UCI')],cntr_disp,greenCI,0.5,4);
    p(3) = shadeplot(realDat.(imgSet{ii}).(region{4}).('all'),[dispDat.(imgSet{ii}).(region{4}).('LCI'); ...
        dispDat.(imgSet{ii}).(region{4}).('UCI')],cntr_disp,redCI,0.5,4);

    set(gca,'fontsize',25,'plotboxaspectratio',[1 1 1],'ylim',[0 5],'xtick',[-2:1:2]);
    legend(p,{'V1','V2','MT'});
    xlabel('Horizontal disparity (\circ)');
    ylabel('Probability density');
    title([imgSet{ii}]);


    fCirc{ii} = figure;
    fCirc{ii}.Position = [1000 1000 650 600];
    hold on;

    p(1) = shadeplot(realDat.(imgSet{ii}).(region{1}).('all'),[dispDat.(imgSet{ii}).(region{1}).('LCI'); ...
        dispDat.(imgSet{ii}).(region{1}).('UCI')],cntr_disp,[0 0 0],0.5,4);

    set(gca,'fontsize',25,'plotboxaspectratio',[1 1 1],'ylim',[0 3.5],'xtick',[-2:1:2]);
    legend(p(1),{'circ'});
    xlabel('Horizontal disparity (\circ)');
    ylabel('Probability density');


end


%% Save figs

saveas(f{1},'./plots/DisparityStats/sando.svg');
saveas(f{2},'./plots/DisparityStats/walking.svg');

saveas(fCirc{1},'./plots/DisparityStats/circSando.svg');
saveas(fCirc{2},'./plots/DisparityStats/circWalking.svg');

