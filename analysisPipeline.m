% This script provides a walk through of the analysis pipeline for the paper. 
% All functions/scripts should be run within the
% top level working directory. After making sure the required datasets are present (see README), 
% run these two commands before running any of the scripts/functions below
% to be sure everything you need is in your path:

addpath(genpath('./sceneStatsAnalysis'));
addpath(genpath('./helper_functions'));

%-------------------------------------------------------------%
% This script reads in the spiking data, processes the data, and fits tuning
% curves. It's slow (a few hours on an M1 Macbook Pro).
% The fits used for the paper are saved as: analysisFiles/physio/fittingResults[area]_final.mat
% If it's re-run, the outputs will be fittingResults[area].mat

main_FitTuningCurves;

% Calls following custom functions:
% - loadDataV1V2.m
% - loadDataMT.m
% - fit1DGabor.m
% - errfun1DGabor.m
% - screen2retDisp.m
% - uniform_sample_in_range.m

% This companion script enables refitting of individual neurons in cases of
% catastophic failure of the main fitting routine:

main_FitTuningCurvesCustom;

%-------------------------------------------------------------%
% This script selects the subset of cells for analysis based on eccentricity and fit quality,
% then plots the RF locations and fitting results of the selected cells
% It saves out a file called fittingResults_processed.mat, which has all the fits plus some
% extra characterizations

main_AssessAndProcessFits;

% creates plots/stats for:
% -Fig 2A
% -Fig 2C
% -Fig 5B
% -Table 2
% -Table 3

% Calls following custom functions:
% - distributionPlot.m (Jonas Dorn, 2017)
% - ColorIt.m
% - subsampleCells.m


%-------------------------------------------------------------%
% This script computes/plots FI across the populations, assuming condition independence.
% Noise models include: fixed Guassian noise (not used in paper),
% Poisson noise, and a per-neuron fitted Gaussian noise model
% It saves out a file called results_populationFI.mat

main_AnalyzeFI;

% creates plots/stats for:
% -Fig 2B
% -Fig 4B (left panel)
% -Table 1

%-------------------------------------------------------------%
% This script explores how various information-limiting correlations between neurons 
% affect the population FI
% It saves out a file called results_populationFI_withCorrelations.mat

main_AnalyzeFI_correlations;

% creates plots for:
% -Fig 4


%-------------------------------------------------------------%
% This script samples disparities from the BORIS image set
% according to sampling distributions derived from neuronal RF center KSDs.
% It is the slowest script, taking a few hours to run in parallel mode on
% an M1 Macbook Pro.
% It saves out a file called dispHistStruct.mat

main_DisparityStats;

% creates plots for:
% -Fig 1B
% -Fig 2D

% Calls following custom functions:
% - dispDistsPlots.m
% - pinky.m (Tristan Ursell, 2012)


%-------------------------------------------------------------%
%This script loads the FI results and compares to the natural disparity histograms
% It saves out a file called results_populationFI_withCorrelations_NDS_comparison.mat

main_ComparePhysioToNDS;

% creates plots/stats for:
% -Fig 3

% Calls following custom function:
% - fitGeneralizedLaplacian


%-------------------------------------------------------------%
% See which parameters of the Gabor fit can best recreate the MT FI
% distribution using the rest of the parameters drawn from V1

main_reparameterizeV1;

% creates plots for:
% -Fig 6
% -Fig 7

% Calls following custom functions:
% - permuteFits.m
% - getJSDiv.m
% - shadeplot.m
% - pinky.m (Tristan Ursell, 2012)


%-------------------------------------------------------------%
% This script makes some additional panels used in Figures

main_AdditionalFigurePanels;

% creates plots for:
% -Fig 1A, C, D
% -Fig 5A
% -Fig 6A


%---------------------------------------------------------------%
% Supplementary analyses:

% Plot the distribution of spike count means/variances between areas 
% to compare "Poisson-ness" for the whole population. This plot is not
% included in the paper, but is included here for reference

supp_FanoFacCheck;

% Plot the population FI for a range of local correlations in which the
% correlation between pairs of neurons is determined by the similarity in
% their preferred disparity

supp_AnalyzeFI_localCorrelations;

% Plot the total population FI for increasing neuronal population sizes to
% confirm that information limiting correlations cause information to
% saturate

supp_AnalyzeFI_infoSaturation;


