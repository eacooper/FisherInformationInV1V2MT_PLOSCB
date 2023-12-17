%% Summarize results of fitting cell populations with Gabors, subsample
% cells using criteria that we want and perform some additional
% calculations used at later stages

clear all; close all;

% Flag = 1 if you want to plot the fit for each cell
plotallfits = 0;

% load symbolic equations for Gabor and first derivate
define_Gabor_tuning_curve;

% define areas
areas = {'V1','V2','MT'};

% define disparity sampling lattice to match disparity statistics
edges_disp = linspace(-2,2,52);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

% Load fitting results
V1 = load('./analysisFiles/physio/fittingResultsV1_final.mat');
V2 = load('./analysisFiles/physio/fittingResultsV2_final.mat');
MT = load('./analysisFiles/physio/fittingResultsMT_final.mat');

% Report number of cells in each population before subsampling
display(['Num V1 units before subsampling: ' num2str(numel(V1.E))]);
display(['Num V2 units before subsampling: ' num2str(numel(V2.E))]);
display(['Num MT units before subsampling: ' num2str(numel(MT.E))]);
display(['Num TOTAL units before subsampling: ' num2str(numel(V1.E)+numel(V2.E)+numel(MT.E))]);


%% Grab RF location data from experiment substructure

% Grab V1 RF location data
for v = 1:length(V1.experiments)
    % positive x is to the animals left in their visual field so we flip to match MT data, positive y is up
    V1.x_pos(v) = -V1.experiments{v}.x_pos;
    V1.y_pos(v) = V1.experiments{v}.y_pos;
end

% Grab V2 RF location data
for v = 1:length(V2.experiments)
    % positive x is to the animals left in their visual field so we flip to match MT data, positive y is up
    V2.x_pos(v) = -V2.experiments{v}.x_pos;
    V2.y_pos(v) = V2.experiments{v}.y_pos;
end

% MT RF location data
for v = 1:length(MT.experiments)
    % positive x is to the animals right in their visual field, positive y is up
    MT.x_pos(v) = MT.experiments{v}.x_pos;
    MT.y_pos(v) = MT.experiments{v}.y_pos;
end


%% Calculate R2 for each fit in each area
% also calculate some other tuning characteristics
xg1     = -2 : 0.01 : 2; % fixed range to evaluate tuning function

for a = 1:length(areas)

    switch areas{a}

        case 'V1';  area = V1;
        case 'V2';  area = V2;
        case 'MT';  area = MT;
        otherwise;  error('invalid area')

    end

    for n = 1:length(area.experiments)

        % grab mean number of repeats per stimulus
        area.mean_repeats(n) = area.experiments{n}.mean_repeats;

        % mean response by disparity
        disps = area.experiments{n}.dat(:,1); % tested disparities
        resps = area.experiments{n}.dat(:,2); % measured responses;

        % calculate fitted spike rates for each tested disparity
        TF = gabor_tuning_curve(disps,area.P(n,:));

        SSR = sum((resps - TF).^2);
        TSS = sum((resps - mean(resps)).^2);

        % R2
        area.r2(n)   = 1 - (SSR/TSS);

        % calculate fitted spike rates for full lattice of disparities
        g1 = gabor_tuning_curve(xg1,area.P(n,:));
        area.allresp(n,:) = g1;

        % max SPS
        area.maxsps(n) = max(g1);

        % preferred disparity
        pref_disp = xg1(g1 == max(g1));
        pref_disp = pref_disp(abs(pref_disp) == min(abs(pref_disp))); % when multiple, take the one with the smallest disparity
        area.pref_disp(n) = pref_disp;

    end

    % fix phase wraparound
    area.P(area.P(:,6)<-pi,6) = area.P(area.P(:,6)<-pi,6) + 2*pi;
    area.P(area.P(:,6)>pi,6)  = area.P(area.P(:,6)>pi,6) - 2*pi;

    % store results
    switch areas{a}

        case 'V1';  V1 = area;
        case 'V2';  V2 = area;
        case 'MT';  MT = area;
        otherwise;  error('invalid area')

    end

end

%% Subsampling:
% -Only keep neurons with R2 > 0.75
% -Crop everything to 10 deg and resample MT to same elevation/abs azimuth as V1/V2
% remove V1/V2 cells where the fixation point was not the center of the
% screen

% Grab indices of cells in each area that meet thess criteria
[subInds] = subsampleCells(V1,V2,MT);

for a = 1:length(areas)

    switch areas{a}

        case 'V1';  area = V1; these_subInds = subInds.V1;
        case 'V2';  area = V2; these_subInds = subInds.V2;
        case 'MT';  area = MT; these_subInds = subInds.MT;
        otherwise;  error('invalid area')

    end

    % keep full datasets for reference
    area_orig = area;

    % now take just the subset of data we want
    area.x_pos   = area.x_pos(these_subInds);
    area.y_pos   = area.y_pos(these_subInds);
    area.allresp = area.allresp(these_subInds,:);
    area.P       = area.P(these_subInds,:);
    area.maxsps  = area.maxsps(these_subInds);
    area.pref_disp = area.pref_disp(these_subInds);
    area.mean_repeats   = area.mean_repeats(these_subInds);
    area.r2   = area.r2(these_subInds);
    area.E = area.E(these_subInds);
    area.S = area.S(these_subInds);
    area.X = area.X(these_subInds);

    area.experiments = [];
    cnt = 1;
    for n = 1:numel(these_subInds)
        if these_subInds(n) == 1
            area.experiments{cnt} = area_orig.experiments{n};
            cnt = cnt + 1;
        end
    end

    switch areas{a}

        case 'V1';  V1 = area; V1orig = area_orig;
        case 'V2';  V2 = area; V2orig = area_orig;
        case 'MT';  MT = area; MTorig = area_orig;
        otherwise;  error('invalid area')

    end

end

%% for each neuron, empirically determine the relationship between
% the mean and variance of the spike count
V1 = fitNeuronalVariances('V1',V1,plotallfits);
V2 = fitNeuronalVariances('V2',V2,plotallfits);
MT = fitNeuronalVariances('MT',MT,plotallfits);
 

%% Plot receptive field locations from the current sample of cells from each area
f1 = figure;
f1.Position = [100 100 650 600];
hold on;
title('Receptive field locations');
scatter(V1.x_pos,V1.y_pos,[],ColorIt('b'),'filled');
scatter(V2.x_pos,V2.y_pos,[],ColorIt('g'),'filled');
scatter(MT.x_pos,MT.y_pos,[],ColorIt('r'),'filled');
plot([-12 12],[0 0],'k:');
plot([0 0],[-12 12],'k:');
axis([-12 12 -12 12]); axis equal tight;
legend('V1','V2','MT');
xlabel('horizontal eccentricity (deg)');
ylabel('vertical eccentricity (deg)');
box on;
set(gca,'fontsize',20);

saveas(f1,'./plots/RFLocations/RFlocations.svg');

%% Perform omnibus + posthoc tests on distributions of parameter fits

% Define area IDs and pairs of areas
pairs  = [1 2; 1 3; 2 3];
pairID = {'V1-V2','V1-MT','V2-MT'};
areaID = {'V1','V2','MT'};
parID  = {'Response offset','Amplitude','Envelope mean','Envelope std.','Frequency','Phase',...
    'RF center eccentricity','Pref. Disparity'};

% Get number of cells in each area
cellCounts = [numel(V1.maxsps) numel(V2.maxsps) numel(MT.maxsps)];

% Make arrays defining indices of concatenated parameter vector
groupVec  = repelem(areaID,cellCounts);
groupInds = repelem([1:3],cellCounts);

%kwTestp = nan(6,1);
%rsTestp = nan(6,3);

% Loop over Gabor parameters
for ii  = 1:6

    % Omnibus non-parametric test (Kruskal-Wallis test)
    % - Concatenate all parameters into a single vector
    if (ii == 3) || (ii == 6)
        % - Since the main difference that affects our hypothesis here is
        % - kurtosis, take the abs for the signed variables
        values{ii}      = abs([V1.P(:,ii)' V2.P(:,ii)' MT.P(:,ii)']);
    else
        values{ii}      = [V1.P(:,ii)' V2.P(:,ii)' MT.P(:,ii)'];
    end

    [kwTestp(ii),thisTable] = kruskalwallis(values{ii},groupVec,'off');
    kwTestStats{ii} = [thisTable{2,3} thisTable{2,5}];

end


% Loop over pairs to do pairwise follow ups
for ii  = 1:6

        % Post-hoc test pairs
    if kwTestp(ii) < 0.05

        % Loop over pairs
        for jj = 1:3

            % Make comparison groups by indexing into vector according to
            % group indices in pairs mat
            [rsTestp(ii,jj),~,thisStruct] = ranksum(values{ii}(groupInds==pairs(jj,1)),values{ii}(groupInds==pairs(jj,2)));
            rsTestStats{ii,jj} = thisStruct.zval;

            % effect size measured as r per Fritz, Morris & Richler 2012 (https://doi.org/10.1037/a0024338)
            rsTestEffectSize{ii,jj} = abs(thisStruct.zval / sqrt(( numel(values{ii}(groupInds==pairs(jj,1))) + numel(values{ii}(groupInds==pairs(jj,2)))) ));

        end

    else

        rsTestp(ii,1) = NaN; rsTestp(ii,2) = NaN; rsTestp(ii,3) = NaN;
        rsTestStats{ii,1} = NaN; rsTestStats{ii,2} = NaN; rsTestStats{ii,3} = NaN;
        rsTestEffectSize{ii,1} = NaN; rsTestEffectSize{ii,2} = NaN; rsTestEffectSize{ii,3} = NaN;

    end

end


% put stats into tables to save and print to workspace
% Define the values for each column
parameters = {'r0', 'A', 'mu', 'sigma', 'f', 'phi'};

% Omnibus Krusal_wallis
for p = 1:6
    chi_sq(p) = round(kwTestStats{p}(2),2);
    pval(p) = kwTestp(p);
end

% Create the table
statsOmni = table(parameters', chi_sq', pval','VariableNames', {'parameter', 'chi_sq', 'p'});

% Display the table
disp('Omnibus Krusal-Wallis test')
disp(statsOmni);

% follow up pairwise
for p = 1:6

    if ~isempty(rsTestStats{p,1})
        V1_V2_Z(p) = round(rsTestStats{p,1},2);
        V1_MT_Z(p) = round(rsTestStats{p,2},2);
        V2_MT_Z(p) = round(rsTestStats{p,3},2);

        V1_V2_r(p) = round(rsTestEffectSize{p,1},2);
        V1_MT_r(p) = round(rsTestEffectSize{p,2},2);
        V2_MT_r(p) = round(rsTestEffectSize{p,3},2);

        V1_V2_pval(p) = rsTestp(p,1);
        V1_MT_pval(p) = rsTestp(p,2);
        V2_MT_pval(p) = rsTestp(p,3);
    else

        V1_V2_Z(p) = NaN; V1_MT_Z(p) = NaN; V2_MT_Z(p) = NaN;

        V1_V2_r(p) = NaN; V1_MT_r(p) = NaN; V2_MT_r(p) = NaN;

        V1_V2_pval(p) = NaN; V1_MT_pval(p) = NaN; V2_MT_pval(p) = NaN;
    end
end

% Create the table
statsPairs = table(parameters', V1_V2_Z', V1_V2_r', V1_V2_pval',...
    V1_MT_Z', V1_MT_r', V1_MT_pval',...
    V2_MT_Z', V2_MT_r', V2_MT_pval',...
    'VariableNames', {'parameter', 'V1_V2_Z', 'V1_V2_r', 'V1_V2_pval', 'V1_MT_Z', 'V1_MT_r', 'V1_MT_pval', 'V2_MT_Z', 'V2_MT_r', 'V2_MT_pval'});

disp('Pairwise Wilcoxon tests');
disp(statsPairs);


%% Plots

% COMPARE PARAMETERS ACROSS REGIONS
% ----------------------- %
% these are the raw fitted parameters, which can be a bit confusing because
% they interact to determine the actual tuning curve shape (e.g., a
% wide/shallow envelope + high amplitude can produce the same max spike
% rate as a narrower envelope and a lower amplitude)
f3 = figure;
f3.Position = [100 100 1500 1000];
hold on;
param_list = {'Response offset','Amplitude','Envelope mean','Envelope std.','Frequency','Phase'};
xlabs = {'Spikes/s','Spikes/s','Disparity (\circ)','Disparity (\circ)','cyc/deg','Radians'};

fTix = [0 1 2 3 4];
for ii = 1:numel(fTix)
    fTixLab{ii} = num2str(fTix(ii));
end

phTix    = [-pi -pi/2 0 pi/2 pi];
phTixLab = {'-\pi' '-\pi/2' '0' '\pi/2' '\pi'};

rtix   = [0.25 10 100 500];
atix   = [2 10 100 500];
envTix = [0.05 0.25 1 4];
for ii = 1:numel(rtix)
    rTixLab{ii}   = num2str(rtix(ii));
    aTixLab{ii}   = num2str(atix(ii));
    envTixLab{ii} = num2str(envTix(ii));
end

for p = 1:6

    subplot(2,3,p);
    hold on;
    title(param_list{p});

    if p == 1 || p == 2 || p == 4
        distributionPlot([log(V1.P(:,p))],'xValues',3,'color',ColorIt('b'),'histOpt',1,...
            'showMM',6,'xNames',{'V1'},'xyOri','flipped');
        distributionPlot([log(V2.P(:,p))],'xValues',2,'color',ColorIt('g'),'histOpt',1,...
            'showMM',6,'xNames',{'V2'},'xyOri','flipped');
        distributionPlot([log(MT.P(:,p))],'xValues',1,'color',ColorIt('r'),'histOpt',1,...
            'showMM',6,'xNames',{'MT'},'xyOri','flipped');

    elseif p == 5
        distributionPlot([(V1.P(:,p))],'xValues',3,'color',ColorIt('b'),'histOpt',0,...
            'divFactor',0:.25:4,'showMM',6,'xNames',{'V1'},'xyOri','flipped')
        distributionPlot([(V2.P(:,p))],'xValues',2,'color',ColorIt('g'),'histOpt',0,...
            'divFactor',0:.25:4,'showMM',6,'xNames',{'V2'},'xyOri','flipped');
        distributionPlot([(MT.P(:,p))],'xValues',1,'color',ColorIt('r'),'histOpt',0,...
            'divFactor',0:.25:4,'showMM',6,'xNames',{'MT'},'xyOri','flipped');

    else
        distributionPlot([V1.P(:,p)],'xValues',3,'color',ColorIt('b'),'histOpt',1,...
            'showMM',6,'xNames',{'V1'},'xyOri','flipped');
        distributionPlot([V2.P(:,p)],'xValues',2,'color',ColorIt('g'),'histOpt',1,...
            'showMM',6,'xNames',{'V2'},'xyOri','flipped');
        distributionPlot([MT.P(:,p)],'xValues',1,'color',ColorIt('r'),'histOpt',1,...
            'showMM',6,'xNames',{'MT'},'xyOri','flipped');

    end

    if p == 1
        set(gca,'xtick',log(rtix),'xticklabel',rTixLab,'xlim',log([0.25 500]));
    elseif p == 2
        set(gca,'xtick',log(atix),'xticklabel',aTixLab,'xlim',log([2 600]));
    elseif p == 4
        set(gca,'xtick',log(envTix),'xticklabel',envTixLab,'xlim',log([0.025 8]));
    elseif p == 5
        set(gca,'xtick',fTix,'xticklabel',fTixLab,'xlim',[0 4]);
    elseif p == 6
        set(gca,'xtick',phTix,'xticklabel',phTixLab,'xlim',[-pi pi]);
    end

    axis equal square; box on;
    set(gca,'ytick',[1 2 3],'YTickLabel',{'MT','V2','V1'},'fontsize',20);
    xlabel(xlabs{p});

end

% Save for Fig 4
saveas(f3,'./plots/AssessFits/Fit_Params.svg');

%% PLOT EXAMPLE FITS
% ----------------------- %
for a = 1:length(areas)

    switch areas{a}

        case 'V1';  area = V1; n = 2; yrange = [30 65];
        case 'V2';  area = V2; n = 2; yrange = [30 80];
        case 'MT';  area = MT; n = 13; yrange = [15 60];
        otherwise;  error('invalid area')

    end

    % plot data and fit
    f5{a} = figure; hold on;
    scatter(area.experiments{n}.dat(:,1),area.experiments{n}.dat(:,2),'k','filled');
    plot(xg1,area.allresp(n,:),'b-','linewidth',2);
    axis square; xticks([]); yticks([]); box on; ylim(yrange); xlim([-2 2]);

    % plot FI using independent Poisson model
    f6{a} = figure; hold on;

    g1 = gabor_tuning_curve(xg1,area.P(n,:));
    dg1 = first_deriv_of_GTC(xg1,area.P(n,:));
    this_FI = dg1.^2 ./ g1;
    plot(xg1,this_FI,'k-','linewidth',2);

    axis square; xticks([]); yticks([]); box on;

    saveas(f5{a},['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '.png']);
    saveas(f5{a},['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '.eps'],'epsc');

    saveas(f6{a},['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '_FI.png']);
    saveas(f6{a},['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '_FI.eps'],'epsc');

end


% PLOT ALL FITS
% ----------------------- %
if(plotallfits)
    for a = 1:length(areas)

        switch areas{a}

            case 'V1';  area = V1;
            case 'V2';  area = V2;
            case 'MT';  area = MT;
            otherwise;  error('invalid area')

        end

        % counters
        pcnt = 1;
        fcnt = 1;
        for n = 1:length(area.experiments)

            if mod(n,66) == 1
                if fcnt > 1
                    saveas(gcf,['./plots/AssessFits/AllFits/Fit_' areas{a} '_' num2str(fcnt-1) '.png']);
                    close gcf;
                end
                figure; hold on;
                setupfig(18,10,10);
                pcnt = 1;
                fcnt = fcnt + 1;
            end

            subplot(6,11,pcnt); hold on; title([ num2str(n) ' R2 = ' num2str(area.r2(n),3)]);

            % data
            scatter(area.experiments{n}.dat(:,1),area.experiments{n}.dat(:,2),'k','filled');

            % our fit
            plot(xg1,area.allresp(n,:),'b-');
            box on;

            xlabel('stimulus disparity (deg)'); ylabel('spike rate mean (sps)');

            pcnt = pcnt + 1;

        end
        saveas(gcf,['./plots/AssessFits/AllFits/Fit_' areas{a} '_' num2str(fcnt-1) '.png']);

    end

end

% save new data for next analysis stage
save('./analysisFiles/physio/fittingResults_processed.mat',...
    'areas','V1','V2','MT','statsOmni','statsPairs','V1orig','V2orig','MTorig','subInds','cntr_disp','edges_disp');
