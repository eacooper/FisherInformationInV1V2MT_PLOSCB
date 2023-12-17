% Assess whether we can account for FI differences between V1 and MT based
% on substituting a single parameter distribution from one to the other

clear all; close all;

% Set the random seed to a specific value so that bootstrapping results are reproducible
rng(724);

% number of times to resample parameters
numResamps = 100;

% load symbolic equations for Gabor and first derivate
define_Gabor_tuning_curve;

% Load in neuronal data
%load('./analysisFiles/physio/results_populationFI.mat');
load('./analysisFiles/results_populationFI_withCorrelations_NDS_comparison.mat');

%% Resample V1 parameter fits from the MT distribution + downsample

% Define parameter names
titles = {'pedestal','amplitude','env. mean','env. width','frequency','phase'};

% for each resampling
for rs = 1:numResamps

    if mod(rs,10) == 0
        disp(['Running resampling run: ',num2str(rs),'/',num2str(numResamps)]);
    end

    % for each Gabor parameter, reassign a value to each V1 cell that is
    % randomly drawn from the parameter probability distribution in MT
    V1_off   = permuteFits(V1,V2,MT,'V1','MT','offset');
    V1_amp   = permuteFits(V1,V2,MT,'V1','MT','amplitude');
    V1_envM  = permuteFits(V1,V2,MT,'V1','MT','envMean');
    V1_envW  = permuteFits(V1,V2,MT,'V1','MT','envWidth');
    V1_freq  = permuteFits(V1,V2,MT,'V1','MT','frequency');
    V1_phase = permuteFits(V1,V2,MT,'V1','MT','phase');

    % perform the same resampling, but draw from the V1 distribution itself
    % as a control
    V1_offCntl   = permuteFits(V1,V2,MT,'V1','V1','offset');
    V1_ampCntl   = permuteFits(V1,V2,MT,'V1','V1','amplitude');
    V1_envMCntl  = permuteFits(V1,V2,MT,'V1','V1','envMean');
    V1_envWCntl  = permuteFits(V1,V2,MT,'V1','V1','envWidth');
    V1_freqCntl  = permuteFits(V1,V2,MT,'V1','V1','frequency');
    V1_phaseCntl = permuteFits(V1,V2,MT,'V1','V1','phase');

    % store the resulting V1 parameters all together
    V1_reparams     = {V1_off.P, V1_amp.P, V1_envM.P, V1_envW.P, V1_freq.P, V1_phase.P};
    V1_reparamsCntl = {V1_offCntl.P, V1_ampCntl.P, V1_envMCntl.P, V1_envWCntl.P, V1_freqCntl.P, V1_phaseCntl.P};

    % Recalulate individual cell FIs under reparameterization

    % For each parameter
    for ii = 1:6

        numCells = size(V1_reparams{ii},1);

        % For each cell
        for jj = 1:numCells

            % grab the resampled and resampled control parameters for this cell
            p  = V1_reparams{ii}(jj,:);
            pC = V1_reparamsCntl{ii}(jj,:);

            % Calculate Fisher information for both
            fi  = first_deriv_of_GTC(cntr_disp,p).^2 ./ gabor_tuning_curve(cntr_disp,p);
            fiC = first_deriv_of_GTC(cntr_disp,pC).^2 ./ gabor_tuning_curve(cntr_disp,pC);

            % Count and remove FI values less than 0, which occur when we
            % got a negative value in the tuning curve
            lt0Cnt_fi  = sum(fi<0);
            lt0Cnt_fiC = sum(fiC<0);

            fi(fi<0)   = 0;
            fiC(fiC<0) = 0;

            % Collect FI distribution into cell array
            V1_reparamsFI{ii,rs}(jj,:)  = fi;
            V1_reparamsFIC{ii,rs}(jj,:) = fiC;

            % also store the number of negative FI values for investigation
            V1_rplt0{ii,rs}(jj)  = lt0Cnt_fi;
            V1_rplt0C{ii,rs}(jj) = lt0Cnt_fiC;

        end

    end
end



%% Characterize the percentage of negative/zero FI values we ended up with 
% to make sure it's not a large amount of the data

res = numel(cntr_disp); % sampling resolution

% for all cells/disparities
perclt0Elems = 100*cellfun(@(x) sum(x),V1_rplt0)/(numCells*(res)); % number of <0 elements for each bootstrap/param
perclt0Cells = 100*cellfun(@(x) sum(x~=0),V1_rplt0)/numCells;      % number of cells with <0 elements for each

perclt0ElemsC = 100*cellfun(@(x) sum(x),V1_rplt0C)/(numCells*(res)); % number of <0 elements for each bootstrap/param
perclt0CellsC = 100*cellfun(@(x) sum(x~=0),V1_rplt0C)/numCells;      % number of cells with <0 elements for each

% averaged for each parameter
meanlt0Elems = mean(perclt0Elems,2);
meanlt0Cells = mean(perclt0Cells,2);
meanlt0ElemsC = mean(perclt0ElemsC,2);
meanlt0CellsC = mean(perclt0CellsC,2);

disp(['Average percent of disparities with FI < 0 in resampled pop. = ' num2str(mean(meanlt0Elems))]);
disp(['Average percent of cells with some FI < 0 in resampled pop. = ' num2str(mean(meanlt0Cells))]);

disp(['Average percent of disparities with FI < 0 in resampled control pop. = ' num2str(mean(meanlt0ElemsC))]);
disp(['Average percent of cells with some FI < 0 in resampled control pop. = ' num2str(mean(meanlt0CellsC))]);


%% Calculate population FI for each reparameterized V1 distribution (MT params and control)

% Get unnormalized population FI from original two populations
mt         = sum(MT.FI_poiss);
v1         = sum(V1.FI_poiss);

% Get AUC for pop FI + normalize by cell count
mtUN       = sum(mt)/size(MT.FI_poiss,1);
v1UN       = sum(v1)/size(V1.FI_poiss,1);

% Normalize pop FI only by cell count
mtcellNorm = mt / size(MT.FI_poiss,1);
v1cellNorm = v1 / size(V1.FI_poiss,1);

% Normalize pop FI by AUC
mt         = mt / sum(mt);
v1         = v1 / sum(v1);

% Initialize matrices for Jenson-Shannon divergences
JSD         = nan(numResamps,6);
JSDcellNorm = nan(numResamps,6);
JSDC        = nan(numResamps,6);

% for each parameter
for ii = 1:6

    % Initialize matrices for FI 
    theseRSPopFIdist         = nan(numResamps,res);
    theseRSPopFIdistCellNorm = nan(numResamps,res);
    theseCntlRSPopFIdist     = nan(numResamps,res);

    % for each resampling that was run
    for rs = 1:numResamps

        % Collect FI distributions for each cell in this bootstrapped
        % sample from MT distribution 
        theseRSFIdists     = V1_reparamsFI{ii,rs};
        theseCntlRSFIdists = V1_reparamsFIC{ii,rs};

        % number of cells
        sampleCellCnt      = size(theseRSFIdists,1);

        % Sum FI over all the cells
        thisRSPopFIdist         = sum(theseRSFIdists);

        % Sum the area under this un-normalized population FI and divide by total cell count
        theseRSTFI(rs,ii)       = sum(thisRSPopFIdist) / sampleCellCnt;

        % Normalize the population FI only by number of cells in sample
        thisRSPopFIdistCellNorm = thisRSPopFIdist / sampleCellCnt;

        % Normalize the population FI so AUC sums to 1
        thisRSPopFIdist         = thisRSPopFIdist / sum(thisRSPopFIdist);

        % Repeat above steps for the bootstraps sampled from V1 distribution
        thisCntlRSPopFIdist         = sum(theseCntlRSFIdists);
        theseCntlRSTFI(rs,ii)       = sum(thisCntlRSPopFIdist) / sampleCellCnt;
        thisCntlRSPopFIdistCellNorm = thisCntlRSPopFIdist / sampleCellCnt;
        thisCntlRSPopFIdist         = thisCntlRSPopFIdist / sum(thisCntlRSPopFIdist);

        % Calculate JSD associated with resampling from MT, both normalized
        % by AUC and by cell num
        JSD(rs,ii)         = getJSDiv(squeeze(thisRSPopFIdist),mt);
        JSDcellNorm(rs,ii) = getJSDiv(squeeze(thisRSPopFIdistCellNorm),mtcellNorm);

        % Calculate JSD for resampling from V1, for AUC normalized data
        JSDC(rs,ii)       = getJSDiv(squeeze(thisCntlRSPopFIdist),mt);

        % Collect into matrices for stats
        theseRSPopFIdist(rs,:)         = thisRSPopFIdist;
        theseRSPopFIdistCellNorm(rs,:) = thisRSPopFIdistCellNorm;
        theseCntlRSPopFIdist(rs,:)     = thisCntlRSPopFIdist;

    end

    % Loop over elements in disparity support and obtain 25th, 50th, and
    % 75th quantiles for plotting
    for jj = 1:(res)

        % Find median FI for this support element and index into KSD support
        % (Sampling from MT)
        v1_med(ii,jj) = median(squeeze(theseRSPopFIdist(:,jj)),'omitnan');
        v1_LEB(ii,jj) = prctile(squeeze(theseRSPopFIdist(:,jj)),[25]);
        v1_UEB(ii,jj) = prctile(squeeze(theseRSPopFIdist(:,jj)),[75]);

        v1_medcellNorm(ii,jj) = median(squeeze(theseRSPopFIdistCellNorm(:,jj)),'omitnan');
        v1_LEBcellNorm(ii,jj) = prctile(squeeze(theseRSPopFIdistCellNorm(:,jj)),[25]);
        v1_UEBcellNorm(ii,jj) = prctile(squeeze(theseRSPopFIdistCellNorm(:,jj)),[75]);

        % (Sampling from V1)
        v1_medC(ii,jj) = median(squeeze(theseCntlRSPopFIdist(:,jj)),'omitnan');
        v1_LEBC(ii,jj) = prctile(squeeze(theseCntlRSPopFIdist(:,jj)),[25]);
        v1_UEBC(ii,jj) = prctile(squeeze(theseCntlRSPopFIdist(:,jj)),[75]);

    end

end


%% Calculate metastatistics

% Sampled from MT
medianJSD     = median(JSD,1);
iqrJSD        = prctile(JSD,[25 75],1);
medianTFI     = median(theseRSTFI,1);
iqrTFI        = prctile(theseRSTFI,[25 75],1);

medianJSDcellNorm     = median(JSDcellNorm,1);
iqrJSDcellNorm        = prctile(JSDcellNorm,[25 75],1);

% Sampled from V1
medianJSDC     = median(JSDC,1);
iqrJSDC        = prctile(JSDC,[25 75],1);
medianTFIC     = median(theseCntlRSTFI,1);
iqrTFIC        = prctile(theseCntlRSTFI,[25 75],1);

% report how much lower the JSD was for the resampled mu population
display(['JSD from resampling mu was these factors lower than the other params:' ...
    num2str(medianJSD/medianJSD(3))]);

display(['JSD from resampling mu was these factors lower than the other params in the control:' ...
    num2str(medianJSDC/medianJSD(3))]);


%% Plot each reparameterized V1 distribution against MT
f1 = figure;
f1.Position = [100 100 1500 1000];
pl = [];

for ii = 1:6
    subplot(2,3,ii);
    hold on;

    % Plot shadeplot
    pl(1) = shadeplot(v1_med(ii,:),[v1_LEB(ii,:); v1_UEB(ii,:)],cntr_disp,[0 0 0],0.4,2.5);
    pl(2) = plot(cntr_disp, mt,'color',ColorIt('r'),'linewidth',2.5);
    pl(3) = plot(cntr_disp, v1,'--k','linewidth',2.5);

    if ii == 1
        legend(pl,{'V1','MT','V1 (empir.)'},'location','northwest');
    end

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:1:2,'ylim',[0 0.2],'yscale','lin','fontsize',15);
    xlabel('Disparity');
    ylabel('Normalized Fisher information');
    title(titles{ii},'AUC normalized');

end

f2 = figure;
f2.Position = [1100 100 1500 1000];
pl = [];

for ii = 1:6
    subplot(2,3,ii);
    hold on;

    % Since we want to plot on log scale here due to huge differences in
    % linear scale between pars, replace lower CI zero values with eps
    %v1_LEBcellNorm(v1_LEBcellNorm==0) = eps;

    % Plot shadeplot
    pl(1) = shadeplot(v1_medcellNorm(ii,:),[v1_LEBcellNorm(ii,:); v1_UEBcellNorm(ii,:)],cntr_disp,[0 0 0],0.4,2.5);
    pl(2) = plot(cntr_disp, mtcellNorm,'color',ColorIt('r'),'linewidth',2.5);
    pl(3) = plot(cntr_disp, v1cellNorm,'--k','linewidth',2.5);
    
    if ii == 1
        legend(pl,{'V1','MT','V1 (empir.)'},'location','northwest');
    end

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:1:2,'ylim',[0 2000],'yscale','lin','fontsize',15);
    %set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:1:2,'ylim',[1E-1 1.9E4],'yscale','log','fontsize',15);
    xlabel('Disparity');
    ylabel('Fisher information normalized by cell count');
    title(titles{ii},' Cell count normalized');

end


%% Compare JS divergence between V1 reparameterizations and MT FI distribution

f3 = figure;
f3.Position = [800 100 625 600];
hold on; pl = [];

for ii = 1:6

    swarmchart(ii*ones(size(JSD,1),1),JSD(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.4);

    pl(1) = scatter(ii,medianJSD(ii),100,'k','filled');
    e = errorbar(ii,medianJSD(ii),medianJSD(ii)-iqrJSD(1,ii),iqrJSD(2,ii)-medianJSD(ii),'k','linewidth',2);

end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log','fontsize',15);
xtickangle(45);
ylabel('JSD');
title('V1-MT Population FI Divergence');

%% Compare JS divergence between V1 reparameterizations and V1 Shuffled FI distribution

f4 = figure;
f4.Position = [800 100 625 600];
hold on; pl = [];

for ii = 1:6

    [p(ii),~,thisStruct] = ranksum(JSD(:,ii),JSDC(:,ii));
    statsJSD(ii,1)       = thisStruct.zval;

    swarmchart(ii*ones(size(JSD,1),1)+0.125,JSD(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.4);
    swarmchart(ii*ones(size(JSDC,1),1)-0.125,JSDC(:,ii),50,[.3 .3 .3] + 0.7*ColorIt('o'),'XJitterWidth',0.4);

    pl(1) = scatter(ii+0.125,medianJSD(ii),100,ColorIt('k'),'filled');
    e = errorbar(ii+0.125,medianJSD(ii),medianJSD(ii)-iqrJSD(1,ii),iqrJSD(2,ii)-medianJSD(ii),'k','linewidth',2);
    e.Color = ColorIt('k');

    pl(2) = scatter(ii-0.125,medianJSDC(ii),100,ColorIt('o'),'filled');
    e = errorbar(ii-0.125,medianJSDC(ii),medianJSDC(ii)-iqrJSDC(1,ii),iqrJSDC(2,ii)-medianJSDC(ii),'k','linewidth',2);
e.Color = ColorIt('o');
end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log','fontsize',15);
xtickangle(45);
ylabel('JSD');
title('V1-MT Population FI Divergence');
legend(pl,{'MT','V1 (shuffled)'},'Location','southeast');


%% Plot AUC ("total fisher info") for the true and reparameterized distributions
f5 = figure;
f5.Position = [800 800 650 600];
hold on; pl = [];

for ii = 1:6

    [p2(ii),~,thisStruct] = ranksum(theseRSTFI(:,ii),theseCntlRSTFI(:,ii));
    statsTFI(ii,1)        = thisStruct.zval;

    swarmchart(ii*ones(size(theseRSTFI,1),1)+0.125,theseRSTFI(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.4);
    swarmchart(ii*ones(size(theseCntlRSTFI,1),1)-0.125,theseCntlRSTFI(:,ii),50,[.3 .3 .3] + 0.7*ColorIt('o'),'XJitterWidth',0.4);

    pl(1) = scatter(ii+0.125,medianTFI(ii),100,ColorIt('k'),'filled');
    e = errorbar(ii+0.125,medianTFI(ii),medianTFI(ii)-iqrTFI(1,ii),iqrTFI(2,ii)-medianTFI(ii),'','linewidth',2);
    e.Color = ColorIt('k');

    pl(2) = scatter(ii-0.125,medianTFIC(ii),100,ColorIt('o'),'filled');
    e = errorbar(ii-0.125,medianTFIC(ii),medianTFIC(ii)-iqrTFIC(1,ii),iqrTFIC(2,ii)-medianTFIC(ii),'k','linewidth',2);
    e.Color = ColorIt('o');

    plot([0 7],mtUN*[1 1],'--r','linewidth',2.5);

end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log',...
    'ylim',5*[1e2 1e6],'fontsize',15);
xtickangle(45);
ylabel('AUC');
title('V1-MT Total Population FI');
legend(pl,{'MT','V1 (shuffled)'},'Location','northeast');


%% Save plots
saveas(f1,'./plots/V1MTReparameterization/reparameterizationBootstraps.svg');
saveas(f2,'./plots/V1MTReparameterization/reparameterizationBootstrapsCellNumNorm.svg');
saveas(f3,'./plots/V1MTReparameterization/reparJSD.svg');
saveas(f4,'./plots/V1MTReparameterization/reparJSD_withV1.svg');
saveas(f5,'./plots/V1MTReparameterization/reparTFI.svg');


