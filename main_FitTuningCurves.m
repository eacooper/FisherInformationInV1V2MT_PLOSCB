% This script reads in the spiking data, processes, and fits tuning
% curves The fits used for the paper are stored with
% results[area]_final.mat, if it's re-run the outputs will be results[area].mat
% - areas: flag which cell population sample we want to fit

clear all; close all;

% you can select one or more areas:
% areas = {'V1'};
% areas = {'V2'};
% areas = {'MT'};
areas = {'V1','V2','MT'};

% V1
if any(contains(areas,'V1'))

    % load and preprocess data
    experiments = loadDataV1V2(1);

    % fitting info we'll store
    P     = []; % parameters
    S     = []; % interpolated spike rates used for fitting
    X     = []; % horizontal disparity used for fitting
    E     = []; % fitting error

    % fit each cell that meets inclusion criteria
    for n = 1 : length(experiments)
        dat = experiments{n}.dat; % pre-processed data

        figure(1); clf;
        [P,S,X,E] = fit1DGabor( 'V1', dat, P, S, X, E, n, length(experiments) );
    end

    save( './analysisFiles/physio/fittingResultsV1.mat', 'P', 'S', 'X', 'E', 'experiments' );

end

% V2
if any(contains(areas,'V2'))

    % load and preprocess data
    experiments = loadDataV1V2(2);

    % fitting info we'll store
    P     = []; % parameters
    S     = []; % interpolated spike rates used for fitting
    X     = []; % horizontal disparity used for fitting
    E     = []; % fitting error

    % fit each cell that meets inclusion criteria
    for n = 1 : length(experiments)
        dat = experiments{n}.dat; % pre-processed data

        figure(1); clf;
        [P,S,X,E] = fit1DGabor( 'V2', dat, P, S, X, E, n, length(experiments) );
    end

    save( './analysisFiles/physio/fittingResultsV2.mat', 'P', 'S', 'X', 'E', 'experiments' );

end

% MT
if any(contains(areas,'MT'))

    % load and preprocess data
    experiments = loadDataMT();

    % fitting info we'll store
    P     = []; % parameters
    S     = []; % interpolated spike rates used for fitting
    X     = []; % horizontal disparity used for fitting
    E     = []; % fitting error

    % fit each cell that meets inclusion criteria
    for n = 1 : length(experiments)
        dat = experiments{n}.dat;

        figure(1); clf;
        [P,S,X,E] = fit1DGabor( 'MT', dat, P, S, X, E, n, length(experiments) );
    end

    save( './analysisFiles/physio/fittingResultsMT.mat', 'P', 'S', 'X', 'E', 'experiments' );
end
