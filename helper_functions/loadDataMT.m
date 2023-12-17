function [experiments] = loadDataMT()

% Get listing of raw data files
D = dir('./dataPhysio/MT/DispTunRaw*.mat');
experiments = {}; % initialize

% load meta-data to grab RF eccentricities
T_MT = readtable('./dataPhysio/MT/metadataMT.xlsx');

cnt = 1;

% for each data file
for k = 1 : numel(D)

    % load it
    load(['./dataPhysio/MT/' D(k).name]);

    % set valid indices for disparity trials
    dispInds = ddat.disparity~=-9999 & ddat.disparity~=-99 & ddat.disparity~=99 & ddat.disparity~=98;

    % horizontal disparity levels and counts for all stimulus presentations
    dx = ddat.disparity(dispInds);
    count = ddat.firing_rate(dispInds);

    % calculate mean spike rate and variance for each unique disparity
    [disparities,~,ic] = unique(dx); % get unique disparities tested and indices for these disparities (indices used below for anova)
    
    % initialize data vectors
    mean_count = zeros(size(disparities));
    var_count = zeros(size(disparities));
    repeats = zeros(size(disparities));

    % for each disparity, count number of repeats, mean spike count, and variance
    for d = 1:length(disparities)
        repeats(d) = sum(dx==disparities(d));
        mean_count(d) = mean(count(dx==disparities(d))); 
        var_count(d) = var(count(dx==disparities(d)));

    end

    % store mean repeats across all disparities
    mean_repeats = mean(repeats);

    % compare filename to metadata in order to get RF eccentricity
    fname = D(k).name;
    fname = fliplr(strtok(fliplr(fname),'_'));
    fname = fname(1:end-4);

    % find row in metadata file
    this_row = find(strcmp(T_MT.FIleName,fname));

    % Grab RF location
    % atan(x/d) where x is displacement on screen, d is distance to nodal point from the screen
    % positive x is right in VF, positive y is up in VF
    x_pos  = T_MT.RF_x(this_row); % note this does not account for vergence, see Jenny's notes
    y_pos  = T_MT.RF_y(this_row);

    % apply screen disparity distortion correction to get disparity in
    % units where 0 = horopter rather than 0 = screen (also, corrects
    % for foreshortening on flat screen)
    horiz_ecc_deg = x_pos;
    disparities = screen2retDisp(disparities,horiz_ecc_deg,'MT');

    % to match V1/V2 data, we run an ANOVA and only keep cells with p <
    % 0.01 (threshold Bruce used to get samples for us)
    groups = ic;
    p  = anova1(count,groups,'off');

    if p >= 0.01
        display(['p >= 0.01 - skipping ' num2str(k)]);
        continue
    end

    % also only keep neurons with average repeats 3 or more
    if mean_repeats < 3
        display(['fewer than 3 repeats on average - skipping ' num2str(k)]);
        continue
    end

    % store into data structure
    experiments{cnt}.dat            = [disparities' mean_count var_count repeats];
    experiments{cnt}.fn             = D(k).name;
    experiments{cnt}.mean_repeats   = mean_repeats;
    experiments{cnt}.x_pos = x_pos;
    experiments{cnt}.y_pos = y_pos;

    cnt = cnt + 1;

end



