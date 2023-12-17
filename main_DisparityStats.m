% This script generates KSDs from the V1, V2, and MT neuronal datasets,
% samples from the true and bootstrapped BORIS image sets,

clear all; close all;

% load neuronal data for plotting RF location KSDs
load('./analysisFiles/physio/fittingResults_processed.mat');

% sampling rate
numSamps = 100; % number of times to sample per disparity map
numBoots = 100; % number of times to bootstrap

% Parallel computing on (not recommended to turn off unless testing)
parallelOn = 1;

% define sampling lattice and eccentricities for disparity maps
rangeMax = 10;         % Limit sampling window to +/-10deg around fixation point
dx       = 10/103;     % interval between each pixel in degrees (from BORIS dataset)

% specify the coordinates of pixels in the BORIS dataset degrees
supp1D   = -rangeMax:dx:rangeMax;
[gx,gy]  = meshgrid(supp1D,supp1D);
gxL      = gx(:);
gyL      = gy(:);
ecc      = sqrt(gx.^2 + gy.^2); % eccentricity in deg
suppSz   = size(gx,1);


%% Generate kernel-smoothed densities from the neuronal RF centers + plot,
% using a sampling lattice that matches to the BORIS dataset

% V1 density plot
V1density = ksdensity([V1.x_pos' V1.y_pos'],[gxL gyL]);
V1densityMat = reshape(V1density,[suppSz suppSz]);

% V2 density plot
V2density = ksdensity([V2.x_pos' V2.y_pos'],[gxL gyL]);
V2densityMat = reshape(V2density,[suppSz suppSz]);

% MT density plot
MTdensity = ksdensity([MT.x_pos' MT.y_pos'],[gxL gyL]);
MTdensityMat = reshape(MTdensity,[suppSz suppSz]);

% Circular center 10 deg
CircDensity = sqrt(gx.^2 + gy.^2)<=10;
circDensityMat = CircDensity/sum(CircDensity(:));

% Save these to plug into image stats script
save(['./analysisFiles/visualFieldSampling/V1densityMat_BORIS.mat'],'V1densityMat','V1');
save(['./analysisFiles/visualFieldSampling/V2densityMat_BORIS.mat'],'V2densityMat','V2');
save(['./analysisFiles/visualFieldSampling/MTdensityMat_BORIS.mat'],'MTdensityMat','MT');
save(['./analysisFiles/visualFieldSampling/circDensityMat_BORIS.mat'],'circDensityMat');


% Plot these KSDs
tix  = [-10 -5 0 5 10];
lims = [-10 10];

cmap = parula;
gam  = 255*linspace(0,1,256).^0.25 + 1;
cmap = [interp1(1:256,cmap(:,1),gam)' interp1(1:256,cmap(:,2),gam)' interp1(1:256,cmap(:,3),gam)'];

% for each area
for a = 1:length(areas)

    this_area = areas{a};

    % grab the data structure for this area
    switch this_area
        case 'V1';      area = V1; densityMat = V1densityMat;
        case 'V2';      area = V2; densityMat = V2densityMat;
        case 'MT';      area = MT; densityMat = MTdensityMat;
    end

    f1 = figure; hold on;
    f1.Position = [100 100 760 600];
    imagesc(tix,tix,densityMat);
    axis image; axis xy;

    scatter(area.x_pos,area.y_pos,30,'w','filled','markerfacealpha',0.5);

    set(gca,'PlotBoxAspectRatio',[1 1 1],'xlim',lims,'ylim',lims,'xtick',tix,'ytick',tix,'fontsize',20);
    xlabel('horizontal eccentricity (deg)');
    ylabel('vertical eccentricity (deg)');
    title(this_area);
    cb1 = colorbar;
    cb1.Label.String = 'probability';
    clim([0 0.06]);
    colormap(cmap);

    saveas(f1,['./plots/RFlocations/' this_area '_KSD.svg']);

end

close all;

%% Sample from BORIS image set and do bootstraps

% Grab disparity images from BORIS dataset
load('./dataSceneStats/Making_sandwich_disparity.mat');
sando = horizontal_disparity;
load('./dataSceneStats/Walking_outside_disparity.mat');
walking = horizontal_disparity;

% sampling resolution
res = numel(edges_disp);

% Get stats from V1/V2/MT from original KSDs + V1/V2 limited MT
getDispStats('sando',sando,numSamps,ecc,edges_disp,res,0,V1densityMat,V2densityMat,MTdensityMat,circDensityMat);
getDispStats('walking',walking,numSamps,ecc,edges_disp,res,0,V1densityMat,V2densityMat,MTdensityMat,circDensityMat);

% Setup parallel processing workers
if parallelOn
    numCores = feature('numcores');
    nCores   = numCores;
    maxNumCompThreads(nCores);
    parpool(nCores);

    parfor ii = 1:numBoots
        resampIter = ii;

        getDispStats('sando',sando,numSamps,ecc,edges_disp,res,resampIter,V1densityMat,V2densityMat,MTdensityMat,circDensityMat);
        getDispStats('walking',walking,numSamps,ecc,edges_disp,res,resampIter,V1densityMat,V2densityMat,MTdensityMat,circDensityMat);

    end
else
    for ii = 1:numBoots
        resampIter = ii;

        getDispStats('sando',sando,numSamps,ecc,edges_disp,res,resampIter,V1densityMat,V2densityMat,MTdensityMat,circDensityMat);
        getDispStats('walking',walking,numSamps,ecc,edges_disp,res,resampIter,V1densityMat,V2densityMat,MTdensityMat,circDensityMat);

    end
end


%% Collect image stats

numHistBins = res-1;

% Collect individual bootstrap stat runs into a single structure

% Define two behaviorally-ID'd image sets
imgSet    = {'sando','walking'};

% Define how we defined the KSDs
KSDid     = {'Circ','V1','V2','MT'};

% Define binning for how we calculate CIs
numCIBins = 10;

% Initialize loop variables
dispDat = struct;
realDat = struct;

% Loop over image sets (sando/walking)
for jj = 1:numel(imgSet)

    % Loop over KSDs
    for kk = 1:numel(KSDid)

        % Load in the real stats
        realDatStruc = load(['./analysisFiles/disparityStats/dispHist',KSDid{kk},'_',imgSet{jj},'.mat']);
        fieldName    = fieldnames(realDatStruc);
        thisDat      = getfield(realDatStruc,fieldName{1});

        realDat = setfield(realDat,imgSet{jj},KSDid{kk},'all',thisDat);

        % Loop over bootstraps
        theseRuns = nan(numBoots,numHistBins);

        for ii = 1:numBoots

            % Collect the data from this run's structure field
            datStruc  = load(['./analysisFiles/disparityStats/bootstraps/dispHist',KSDid{kk},'_',imgSet{jj},num2str(ii),'.mat']);
            fieldName = fieldnames(datStruc);
            thisRun   = getfield(datStruc,fieldName{1});

            theseRuns(ii,:) = thisRun;
        end

        % Pack this data into a new structure field
        dispDat = setfield(dispDat,imgSet{jj},KSDid{kk},'all',theseRuns);

        % Get 95% CIs for each of the histogram bins
        alpha = 0.05;

        for ii = 1:numHistBins

            [thisKSD,theseVals] = ksdensity(theseRuns(:,ii));

            thisCumDens = cumsum(thisKSD/sum(thisKSD));

            [~,LCIind] = min(abs(thisCumDens - alpha/2));
            [~,UCIind] = min(abs(thisCumDens - (1 - alpha/2)));

            theseLCI(ii) = theseVals(LCIind);
            theseUCI(ii) = theseVals(UCIind);

        end

        dispDat = setfield(dispDat,imgSet{jj},KSDid{kk},'LCI',theseLCI);
        dispDat = setfield(dispDat,imgSet{jj},KSDid{kk},'UCI',theseUCI);

    end

end

% Save the collected stats data
save('./analysisFiles/disparityStats/dispHistStruct.mat','dispDat','realDat','cntr_disp');

%% Plot disparity probability distributions used for main figs in paper and controls

% Load in disparity histograms based on each KSD (keep this here in case
% you don't want to rerun the first chunk of this script since it's slow
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


