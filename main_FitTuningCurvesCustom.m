% Handle a couple of edge cases for tuning curve fitting
% To refit an individual cell, the area and the cell number should be uncommented

clear all; close all;

%area = 'V1';  % max fitted spike rate exceeds measured max by a factor of
%10, better fits acheived penalizing fits where the max fitted rate is more than
%2x the observed max

%n = 70;
%n = 74;
%n = 195;

%area = 'V2'; % max fitted spike rate exceeds measured max by a factor of
%10, better fits acheived penalizing fits where the max fitted rate is more than
%2x the observed max

%n = 485; 
%n = 513; 


%area = 'MT'; % gaps between disparity samples are too broad so we got fits
% with aliasing (high frequency that goes up and down through all points.
% Re-running with max freq as 3 cpd

%n = 15;
%n = 19;
%n = 52;
%n = 116;
%n = 120;
%n = 157;
%n = 164;
%n = 165;
%n = 171;
%n = 181;
%n = 185;
%n = 210;
%n = 232;
%n = 240;

%n = 250;
%n = 293;
%n = 294;
%n = 325;
%n = 373;

% Re-running with max freq as 1 cpd
%n = 245;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flag indicating whether to correct the stimulus disparities to account
% for planar screen (horopter deviation and foreshortening)

correct_screen_disparity = 1;

% load and preprocess data
if strcmp(area,'V1')
    experiments = loadDataV1V2(1);
elseif strcmp(area,'V2')
    experiments = loadDataV1V2(2);
elseif strcmp(area,'MT')
    experiments = loadDataMT();
end

% fitting info we'll store
P     = []; % parameters
S     = []; % interpolated spike rates used for fitting
X     = []; % horizontal disparity used for fitting
E     = []; % fitting error

dat = experiments{n}.dat; % pre-processed data for this cell

figure(1); clf;

% do the fitting for just this cell
[P,S,X,E] = fit1DGabor( area, dat, P, S, X, E, 1, 1 );

% get cell index
Pind = P; Sind = S; Xind = X; Eind = E;

% load in full data
if strcmp(area,'V1')
    load('./analysisFiles/physio/fittingResultsV1.mat');
elseif strcmp(area,'V2')
    load('./analysisFiles/physio/fittingResultsV2.mat');
elseif strcmp(area,'MT')
    load('./analysisFiles/physio/fittingResultsMT.mat');
end

% replace this cell with refit
P(n,:)   = Pind;
S(n).s   = Sind.s;
X(n).x   = Xind.x;
E(n)     = Eind;

% save results
if strcmp(area,'V1')
    save( './analysisFiles/physio/fittingResultsV1.mat', 'P', 'S', 'X', 'E', 'experiments' );
elseif strcmp(area,'V2')
    save( './analysisFiles/physio/fittingResultsV2.mat', 'P', 'S', 'X', 'E', 'experiments' );
elseif strcmp(area,'MT')
    save( './analysisFiles/physio/fittingResultsMT.mat', 'P', 'S', 'X', 'E', 'experiments' );
end
