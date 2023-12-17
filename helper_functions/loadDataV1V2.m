function [experiments] = loadDataV1V2(area)

% Load in raw data
load('./dataPhysio/V1V2/AllDTData.mat');

experiments = {};
cnt = 1;

% for each entry
for x = 1:length(DTdata)

    % if this entry belongs to the right area
    if DTdata(x).visualarea == area

        % load data
        disparities = DTdata(x).dx;
        mean_count  = DTdata(x).resp;
        var_count   = cellfun(@(z) var(z),DTdata(x).counts);

        % some data has more than 2 dimensions
        var_count = squeeze(var_count);

        % count number of repeats per disparity and mean repeats across all disparities
        counts = DTdata(x).counts;
        repeats = [];
        for c = 1:length(counts)
            repeats(c) = length(DTdata(x).counts{c});
        end
        mean_repeats = mean(repeats);

        % grab RF location
        % atan(x/d) where x is displacement on screen, d is distance to nodal point from the screen
        % positive x is to the animals left in their visual field, positive y is up
        x_pos = DTdata(x).rf(1);
        y_pos = DTdata(x).rf(2);

        % Note: These data were already sampled from a larger dataset by running an ANOVA 
        % and only keeping cells with p < 0.01

        % remove non disparity conditions
        remove      = disparities < -500;
        disparities = disparities(~remove);
        mean_count  = mean_count(~remove)';
        var_count   = var_count(~remove)';
        repeats     = repeats(~remove);

        % ignore cases where there are repeated mean responses given for the same disparity
        if numel(unique(disparities)) < numel(disparities)
            display(['unresolved repeats - skipping ' num2str(x)]);
            continue
        end

        % apply screen disparity distortion correction to get disparity in
        % units where 0 = horopter rather than 0 = screen (also, corrects
        % for foreshortening on flat screen)
        % NOTE: this transformation is only valid for cells where the
        % fixation point (.fp) is 0,0 or NaN,NaN, because these were
        % centered on the screen. In practice, we run this calculation on
        % everything and filter based on fp later on.
        horiz_ecc_deg = x_pos;
        disparities = screen2retDisp(disparities,horiz_ecc_deg,['V' num2str(area)]);
        disparities = disparities';

        % ignore cases where orientation of disparity is not horizontal
        if ~isnan(DTdata(x).or) && abs(DTdata(x).or) ~= 90
            display(['disparity not horizontal - skipping ' num2str(x)]);
            continue;
        end

        % when RF orientation is -90 the disparity is still horizontal, but the direction is reversed.  So now negative values are uncrossed.
        if DTdata(x).or == -90
            display(['disparity axis reversed, fixing values ' num2str(x)]);
            disparities = -disparities;
        end

        % also only keep neurons with average repeats 3 or more
        if mean_repeats < 3
            display(['fewer than 3 repeats on average - skipping ' num2str(x)]);
            continue;
        end

        % omit cells fit to a very small range of disparities,
        % because this is hard to compare to MT data
        if range(disparities) < 1
            display(['disparity range too small - skipping ' num2str(x)]);
            continue;
        end


        % convert responses to spikes per second
        % .resp is in units of mean spike count, not divided by duration. Duration was 0.5sec  0.47 or 0.42 sec (different studies).  ]
        % There is a field in the dxresp that records a  number related to the stimulus duration.
        % (This gets complicated, on some trials the graphics card fails to keep up, so it takes 43 video frames to paint 42 actual frames.
        % The duration number in these files depends on the distribution of real durations in the expt.)
        % But I always just count spikes over either 420ms o470, or  500ms, whichever is the nominal (and so shortest) duration.
        % SO if "duration" is <  0.47, its 420ms.  if its > 0.47 and < 0.49, its 470ms (names all begin 'jbeG')  and if its > 0.49, its 500ms (names all begin 'jbeM' or 'lemM')
        if DTdata(x).duration < 0.47
            duration = 0.42;
        elseif DTdata(x).duration >= 0.47 && DTdata(x).duration < 0.49
            duration = .47;
        elseif DTdata(x).duration >= 0.49
            duration = .5;
        else
            error('invalid duration');
        end

        mean_count = mean_count/duration;
        var_count = var_count/duration;

        % store in structure:
        experiments{cnt}.dat = [disparities mean_count' var_count' repeats']; % disparities, responses, and repeats
        experiments{cnt}.mean_repeats = mean_repeats; % mean repeats
        experiments{cnt}.x_pos = x_pos; % RF location
        experiments{cnt}.y_pos = y_pos;
        experiments{cnt}.alldata = DTdata(x); % all other data about this cell

        cnt = cnt + 1;

    end

end


