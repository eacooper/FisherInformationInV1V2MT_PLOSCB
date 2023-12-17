%% This script analyzes (plots FI distribution and RF centers for) the FI distributions
% from each of the cortical areas. We run the basic calculation for Gaussian, Poisson, and Gaussian variable
% noise with no correlations

clear all; close all;

% Set the random seed to a specific value so that bootstrapping results are reproducible
rng(724);

% bootstrapping parameters
iterations  = 100; % Iterations for bootstrapping
sams        = 200; % Sample size for each bootstrapped sample

% Load fitting results that are subsampled to the neurons that we want to analyze
load('analysisFiles/physio/fittingResults_processed.mat');

% load symbolic equations for Gabor and first derivate
define_Gabor_tuning_curve;

% for each area, calculate the values we need to computer FI
for a = 1:length(areas)

    switch areas{a}

        case 'V1';  area = V1;
        case 'V2';  area = V2;
        case 'MT';  area = MT;
        otherwise;  error('invalid area')

    end

    % for each neuron
    for n = 1:length(area.experiments)

        % evaluate tuning curve on disparity stats lattice
        area.mean_responses(n,:) = gabor_tuning_curve(cntr_disp,area.P(n,:));

        % first derivative
        area.first_derivs(n,:) = first_deriv_of_GTC(cntr_disp,area.P(n,:));

        % compute the variance estimate at each disparity given the fitted mean/var line
        area.est_variances(n,:) = area.mean_var_slope(n)*area.mean_responses(n,:) + area.mean_var_intercept(n);

    end

    switch areas{a}

        case 'V1';  V1 = area;
        case 'V2';  V2 = area;
        case 'MT';  MT = area;
        otherwise;  error('invalid area')

    end

end


%% Now we'll just compute the FI for Gaussian, Poisson, and Gaussian variable noisy independent
% populations. These look really similar to each other, which is good!

figGaussUni = figure;
figGaussUni.Position = [100 100 650 600];
figGaussUni.Renderer = 'painter';

figPoiss = figure;
figPoiss.Position = [100 100 650 600];
figPoiss.Renderer = 'painter';

figGaussVar = figure;
figGaussVar.Position = [100 100 650 600];
figGaussVar.Renderer = 'painter';

plot_colors = {'b','g','r'};

for a = 1:length(areas)

    this_area = areas{a};

    % grab the data structure for this area
    switch this_area
        case 'V1';      area = V1;
        case 'V2';      area = V2;
        case 'MT';      area = MT;
    end

    % FI for Gaussian noise assuming uniform variance

    % for each neuron this is just the first derivative squared divided by
    % a constant that reflects the expected variance
    area.FI_gauss_uni_var = mean(area.est_variances(area.est_variances>=0));
    area.FI_gauss_uni = (area.first_derivs.^2)/area.FI_gauss_uni_var;

    % sum up these to get our measure of FI with Gaussian noise and independent neurons
    area.FI_gauss_uni_pop = sum(area.FI_gauss_uni);

    % normalize for each population
    area.FI_gauss_uni_pop_total = sum(area.FI_gauss_uni_pop);
    area.FI_gauss_uni_pop = area.FI_gauss_uni_pop./area.FI_gauss_uni_pop_total;


    % FI for Poisson noise

    % reassign any very small means to the 5th percentile of this area to avoid unstable calculations
    area.mean_responses_clean = area.mean_responses;
    area.mean_responses_clean(area.mean_responses <= quantile(area.mean_responses(:),0.05)) = quantile(area.mean_responses(:),0.05);

    % for each neuron this is the first derivative squared divided by the mean
    area.FI_poiss = area.first_derivs.^2 ./area.mean_responses_clean;

    % sum up to get the population FI with Poisson noise and independent neurons
    % and normalize for each population
    area.FI_poiss_pop = sum(area.FI_poiss);
    area.FI_poiss_pop_total = sum(area.FI_poiss_pop);
    area.FI_poiss_pop = area.FI_poiss_pop./sum(area.FI_poiss_pop_total);


    % FI for Gaussian noise with variance fitted to each cell

    % reassign any very small variances to the 5th percentile of this area
    % to avoid unstable calculations
    area.est_variances_clean = area.est_variances;
    area.est_variances_clean(area.est_variances <= quantile(area.est_variances(:),0.05)) = quantile(area.est_variances(:),0.05);

    % for each neuron this is just the first derivative squared
    area.FI_gauss_var = (area.first_derivs.^2) ./area.est_variances_clean;

    % sum up the squares of these to get our control measure of FI with Gaussian noise and independent neurons
    % normalize for each population
    area.FI_gauss_var_pop = sum(area.FI_gauss_var);
    area.FI_gauss_var_pop_total = sum(area.FI_gauss_var_pop);
    area.FI_gauss_var_pop = area.FI_gauss_var_pop./area.FI_gauss_var_pop_total;

    % add bootstraps
    for n = 1:iterations

        if ~mod(n,50)
            display(['running bootstraps iteration ' num2str(n)]);
        end

        % sample with replacement
        littleFI                    = area.FI_gauss_uni(randsample(size(area.FI_gauss_uni,1),sams,true),:);
        area.FI_gauss_uni_boot(n,:) = sum(littleFI);
        area.FI_gauss_uni_boot(n,:) = area.FI_gauss_uni_boot(n,:) / sum(area.FI_gauss_uni_boot(n,:));

        littleFI                    = area.FI_poiss(randsample(size(area.FI_poiss,1),sams,true),:);
        area.FI_poiss_boot(n,:)     = sum(littleFI);
        area.FI_poiss_boot(n,:)     = area.FI_poiss_boot(n,:) / sum(area.FI_poiss_boot(n,:));

        littleFI                    = area.FI_gauss_var(randsample(size(area.FI_gauss_var,1),sams,true),:);
        area.FI_gauss_var_boot(n,:) = sum(littleFI);
        area.FI_gauss_var_boot(n,:) = area.FI_gauss_var_boot(n,:) / sum(area.FI_gauss_var_boot(n,:));

    end


    % plot it
    figure(figGaussUni); hold on;
    plot(cntr_disp, area.FI_gauss_uni_boot, 'color', 0.5*ones(1,3) + 0.5*ColorIt(plot_colors{a}),'linewidth',0.25);
    plot(cntr_disp, area.FI_gauss_uni_pop,'color',ColorIt(plot_colors{a}),'linewidth',4);

    title('Fixed Gaussian Noise, Independent Neurons')
    axis square; box on; ylim([0 0.23])
    xlabel('horizontal disparity (deg)');
    ylabel('Fisher Information (normalized)');
    set(gca,'fontsize',20);

    figure(figPoiss); hold on;
    plot(cntr_disp, area.FI_poiss_boot, 'color', 0.5*ones(1,3) + 0.5*ColorIt(plot_colors{a}),'linewidth',0.25);
    hP(a) = plot(cntr_disp, area.FI_poiss_pop,'color',ColorIt(plot_colors{a}),'linewidth',4);

    title('Poisson Noise, Independent Neurons')
    axis square; box on; ylim([0 0.23])
    xlabel('horizontal disparity (deg)');
    ylabel('Fisher Information (normalized)');
    set(gca,'fontsize',20);

    figure(figGaussVar); hold on;
    plot(cntr_disp, area.FI_gauss_var_boot, 'color', 0.5*ones(1,3) + 0.5*ColorIt(plot_colors{a}),'linewidth',0.25);
    plot(cntr_disp, area.FI_gauss_var_pop,'color',ColorIt(plot_colors{a}),'linewidth',4);

    title('Fitted Gaussian Noise, Independent Neurons')
    axis square; box on; ylim([0 0.23])
    xlabel('horizontal disparity (deg)');
    ylabel('Fisher Information (normalized)');
    set(gca,'fontsize',20);  


    % store results for each area in a separate matrix
    switch this_area
        case 'V1';      V1 = area;
        case 'V2';      V2 = area;
        case 'MT';      MT = area;
    end

end

figure(figPoiss); legend(hP,'V1','V2','MT');

% plot all cells Poisson FI as heat map
fAll = figure; hold on;
subplot(1,3,1); hold on; title('V1');
imagesc(V1.FI_poiss); colorbar; axis tight;
subplot(1,3,2); hold on; title('V2');
imagesc(V2.FI_poiss); colorbar; axis tight;
subplot(1,3,3); hold on; title('MT');
imagesc(MT.FI_poiss); colorbar; axis tight;

% plot a box and whisker of the slope and intercepts for variance fits
figFits = figure; hold on;
figFits.Position = [100 100 400 200];
subplot(1,2,1); hold on; title('slopes');
hB = boxplot([V1.mean_var_slope V2.mean_var_slope MT.mean_var_slope],...
    [repmat(1,1,numel(V1.mean_var_slope)) repmat(2,1,numel(V2.mean_var_slope)) repmat(3,1,numel(MT.mean_var_slope))],...
    'symbol', '','colors',[ColorIt('b') ; ColorIt('g') ; ColorIt('r')]);
plot([0 4],[1 1],'k-');
ylim([-5 7]); xlim([0.5 3.5]); xticklabels({'V1','V2','MT'})
yticks([-4 -2 0 2 4 6])

subplot(1,2,2); hold on; title('intercepts');
boxplot([V1.mean_var_intercept V2.mean_var_intercept MT.mean_var_intercept],...
    [repmat(1,1,numel(V1.mean_var_slope)) repmat(2,1,numel(V2.mean_var_slope)) repmat(3,1,numel(MT.mean_var_slope))],...
    'symbol', '','colors',[ColorIt('b') ; ColorIt('g') ; ColorIt('r')]);
plot([0 4],[0 0],'k-');
ylim([-210 210]); xlim([0.5 3.5]); xticklabels({'V1','V2','MT'})
yticks([-200 -100 0 100 200])

% also report in a table
slopes = [median(V1.mean_var_slope) quantile(V1.mean_var_slope,.25) quantile(V1.mean_var_slope,.75) ;...
    median(V2.mean_var_slope) quantile(V2.mean_var_slope,.25) quantile(V2.mean_var_slope,.75) ;...
    median(MT.mean_var_slope) quantile(MT.mean_var_slope,.25) quantile(MT.mean_var_slope,.75) ];

intercepts = [median(V1.mean_var_intercept) quantile(V1.mean_var_intercept,.25) quantile(V1.mean_var_intercept,.75) ;...
    median(V2.mean_var_intercept) quantile(V2.mean_var_intercept,.25) quantile(V2.mean_var_intercept,.75) ;...
    median(MT.mean_var_intercept) quantile(MT.mean_var_intercept,.25) quantile(MT.mean_var_intercept,.75) ];

% Create the table
statsMeanVarFits = table(areas', slopes(:,1), slopes(:,2), slopes(:,3),...
    intercepts(:,1), intercepts(:,2), intercepts(:,3),...
    'VariableNames', {'area', 'slope_med', 'slope25', 'slope75', 'intercept_med', 'intercept25', 'intercept75'});

% Display the table
disp('best fit slopes and intercepts for mean/variance relationship')
disp(statsMeanVarFits);

% save these plots
saveas(fAll,'./plots/AnalyzeFI/FI_AllNeurons.png');
saveas(figPoiss,'./plots/AnalyzeFI/FI_IndependentPoisson.svg');
saveas(figGaussUni,'./plots/AnalyzeFI/FI_IndependentGaussianUniform.svg');
saveas(figGaussVar,'./plots/AnalyzeFI/FI_IndependentGaussianFitted.svg');
saveas(figPoiss,'./plots/AnalyzeFI/FI_IndependentPoisson.png');
saveas(figGaussUni,'./plots/AnalyzeFI/FI_IndependentGaussianUniform.png');
saveas(figGaussVar,'./plots/AnalyzeFI/FI_IndependentGaussianFitted.png');

saveas(figFits,'./plots/AnalyzeFI/FI_Variance_Fits.svg');


% Save population FI and bootstraps
save('./analysisFiles/physio/results_populationFI.mat',...
    'areas','V1','V2','MT','statsOmni','statsPairs','statsMeanVarFits','V1orig','V2orig','MTorig','subInds','cntr_disp','edges_disp','iterations','sams');


