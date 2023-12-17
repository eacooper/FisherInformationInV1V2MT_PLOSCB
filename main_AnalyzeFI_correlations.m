%% In the main analysis, we computed Fisher information for each neuron, 
% then summed them up because we assumed
% statistical independence. In this control analysis, we'll test how the
% results might change if the neurons has statistical dependences

clear all; close all;

% Load Gabor fitting results of the neurons that we want to analyze
load('analysisFiles/physio/results_populationFI.mat');

% load symbolic equations for Gabor and first derivate
define_Gabor_tuning_curve;

% flag indicating whether to plot individual information limiting correlation matrices
plot_info_matrices = 1;

% create identity matrix and diagonal covariance matrix for each area
for a = 1:length(areas)

    this_area = areas{a};

    % grab the data structure for this area
    switch this_area
        case 'V1';      area = V1;
        case 'V2';      area = V2;
        case 'MT';      area = MT;
    end

    % how many neurons in this area?
    numNeurons = size(area.allresp,1);

    % make identity matrix
    area.C_identity_corr = eye(numNeurons);

    % for each neuron
    for i = 1:numNeurons
        % for each disparity, convert correlations to variance
        for d = 1:length(cntr_disp)
            area.C_identity(i,i,d) = area.est_variances_clean(i,d);
        end
    end

    switch this_area
        case 'V1';      V1 = area;
        case 'V2';      V2 = area;
        case 'MT';      MT = area;
    end

end


%% Run sanity check for independent Gaussian noise, should match main Gaussian noise calculation

for a = 1:length(areas)

    this_area = areas{a};

    % grab the data structure for this area
    switch this_area
        case 'V1';      area = V1;
        case 'V2';      area = V2;
        case 'MT';      area = MT;
    end

    % sanity check using identity matrix
    for d = 1:length(cntr_disp)
        thisC = area.C_identity(:,:,d);
        invC = inv(thisC);
        area.FI_gauss_identity(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);
    end
    area.FI_gauss_identity_total = sum(area.FI_gauss_identity);
    area.FI_gauss_identity =  area.FI_gauss_identity / area.FI_gauss_identity_total;

    % store results for each area in a separate matrix
    switch this_area
        case 'V1';      V1 = area;
        case 'V2';      V2 = area;
        case 'MT';      MT = area;
    end

end

%% plot the population FIs with Gaussian noise (sanity check)
figure; hold on
plot(cntr_disp, V1.FI_gauss_identity,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_gauss_identity,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_gauss_identity,'color',ColorIt('r'),'linewidth',4);

title('Gaussian Noise, Independent Neurons Sanity Check')
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Fisher Information (normalized)');
set(gca,'fontsize',20); ylim([0 0.2]);
saveas(gcf,'./plots/NoiseCorrelations/FI_IndependentGaussianCheck.png');
saveas(gcf,'./plots/NoiseCorrelations/FI_IndependentGaussianCheck.svg');

%% now add in information limiting correlations, for which the correlation
% matrix is determined by the product of the tuning curve derivatives

% we will iterate over scale factors for each area to select ones that reduce the total
% FI by these amounts:
FI_reduction_ratio_target = [1.25 1.5 2];

% for each area
for a = 1:length(areas)

    this_area = areas{a};

    % grab the data structure for this area and specity custom range for
    % scale factors we will apply
    switch this_area
        case 'V1';      area = V1; area.scale_factors = linspace(0.000001,0.00002,10);
        case 'V2';      area = V2; area.scale_factors = linspace(0.000001,0.00001,10);
        case 'MT';      area = MT; area.scale_factors = linspace(0.00003,0.0002,10);
    end

    for s = 1:length(area.scale_factors)

        % initialize matrix for information limiting correlations
        area.C_info = [];

        for d = 1:length(cntr_disp)

            % construct information limiting correlation matrices, in which
            % pairwise correlations are prop to the product of neuronal derivatives
            % per Moreno-Bore et al. Nat Neuro 2014 pg 1411-1412
            C_info = area.first_derivs(:,d)*area.first_derivs(:,d)';

            C_info = C_info * area.scale_factors(s);

            % add these information limiting correlations to the
            % identity covariance matrix, while keeping the diagonal unchanged
            C_info(area.C_identity_corr(:,:) == 1) = 0;

            C_info = area.C_identity(:,:,d) + C_info;

            % apply information limiting correlations
            invC = inv(C_info);
            FI_gauss_info(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);

        end

        % how did these correlations affect total information?
        FI_gauss_info_total(s) = sum(FI_gauss_info);

    end

    % compute information ratio
    area.FI_reduction_ratio = area.FI_gauss_identity_total./FI_gauss_info_total;

    % store results for each area in a separate matrix
    switch this_area
        case 'V1';      V1 = area;
        case 'V2';      V2 = area;
        case 'MT';      MT = area;
    end

end

% find desired scales
V1.scale_factors_selected = interp1(V1.FI_reduction_ratio,V1.scale_factors,FI_reduction_ratio_target,'linear');
V2.scale_factors_selected = interp1(V2.FI_reduction_ratio,V2.scale_factors,FI_reduction_ratio_target,'linear');
MT.scale_factors_selected = interp1(MT.FI_reduction_ratio,MT.scale_factors,FI_reduction_ratio_target,'linear');

% plot ratios associated with these scale factors
figure; hold on; sgtitle('relationship between scale factor and information reduction')
subplot(1,3,1); hold on; title('V1')
plot(V1.scale_factors, V1.FI_reduction_ratio,'-');
plot(V1.scale_factors_selected, FI_reduction_ratio_target,'o');
xlabel('scale factor'); ylabel('information reduction ratio');

subplot(1,3,2); hold on; title('V2')
plot(V2.scale_factors, V2.FI_reduction_ratio,'-');
plot(V2.scale_factors_selected, FI_reduction_ratio_target,'o');
xlabel('scale factor'); ylabel('information reduction ratio');

subplot(1,3,3); hold on; title('MT')
plot(MT.scale_factors, MT.FI_reduction_ratio,'-');
plot(MT.scale_factors_selected, FI_reduction_ratio_target,'o');
xlabel('scale factor'); ylabel('information reduction ratio');

% sort neurons by their preferred disparity, for plotting correlation matrices
[~,V1I] = sort(V1.pref_disp);
[~,V2I] = sort(V2.pref_disp);
[~,MTI] = sort(MT.pref_disp);

%% apply these correlations
for a = 1:length(areas)

    this_area = areas{a};

    % grab the data structure for this area
    switch this_area
        case 'V1';      area = V1; areaI = V1I;
        case 'V2';      area = V2; areaI = V2I;
        case 'MT';      area = MT; areaI = MTI;
    end

    % specify the scale factor for weak and strong information limiting
    % correlations, which are proportionate to the product of neuronal
    % derivatives. The derivate products end up being large, so these scale
    % factors are quite small
    info_weak_scale     = area.scale_factors_selected(1);
    info_med_scale      = area.scale_factors_selected(2);
    info_strong_scale   = area.scale_factors_selected(3);

    % for each disparity
    for d = 1:length(cntr_disp)

        % construct information limiting correlation matrices, in which
        % pairwise correlations are prop to the product of neuronal derivatives
        % per Moreno-Bore pg 1411-1412
        C_info = area.first_derivs(:,d)*area.first_derivs(:,d)';

        %preserve identity along diagonal
        C_info(area.C_identity_corr(:,:) == 1) = 0;

        C_info_weak = C_info * info_weak_scale;
        C_info_med = C_info * info_med_scale;
        C_info_strong = C_info * info_strong_scale;

        % add these information limiting correlations to the
        % indepdenent covariance matrix, while keeping the diagonal
        % unchanged
        C_info_weak = area.C_identity(:,:,d) + C_info_weak;
        C_info_med = area.C_identity(:,:,d) + C_info_med;
        C_info_strong = area.C_identity(:,:,d) + C_info_strong;

        % store for reference
        area.C_info_weak(:,:,d) = C_info_weak;
        area.C_info_med(:,:,d) = C_info_med;
        area.C_info_strong(:,:,d) = C_info_strong;

        % apply weak information limiting correlations
        invC = inv(C_info_weak);
        area.FI_gauss_info_weak(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);

        % apply medium information limiting correlations
        invC = inv(C_info_med);
        area.FI_gauss_info_med(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);

        % apply strong information limiting correlations
        invC = inv(C_info_strong);
        area.FI_gauss_info_strong(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);

    end

    % weak information limiting correlations
    area.FI_gauss_info_weak_total = sum(area.FI_gauss_info_weak);
    area.FI_gauss_info_weak =  area.FI_gauss_info_weak / area.FI_gauss_info_weak_total;

    % medium information limiting correlations
    area.FI_gauss_info_med_total = sum(area.FI_gauss_info_med);
    area.FI_gauss_info_med =  area.FI_gauss_info_med / area.FI_gauss_info_med_total;

    % strong information limiting correlations
    area.FI_gauss_info_strong_total = sum(area.FI_gauss_info_strong);
    area.FI_gauss_info_strong =  area.FI_gauss_info_strong / area.FI_gauss_info_strong_total;

    % for reference, plot the correlation matrices for a sampling of disparities
    % convert covariances back to correlations by dividing by the diagonal
    % of the identity covariance
    if(plot_info_matrices)
        figCorrs = figure; hold on; sgtitle([this_area ' information limiting matrices']);
        set(gcf,'Position',[10 10 1400 1000]);

        cnt = 1;
        for d = [ 1 13 24 26 28 39 51]
            figure(figCorrs);
            subplot(7,6,cnt); hold on; title(['weak corrs, disparity = ' num2str(cntr_disp(d),2) 'deg']);
            imagesc(area.C_info_weak(areaI,areaI,d) ./ area.est_variances_clean(areaI,d));
            axis image; colorbar;

            subplot(7,6,cnt+1); hold on;
            histogram(area.C_info_weak(:,:,d) ./ area.est_variances_clean(:,d),linspace(-1,1,25));
            xlabel('correlation'); ylabel('log freq'); set(gca,'YScale','log');

            subplot(7,6,cnt+2); hold on; title(['med corrs, disparity = ' num2str(cntr_disp(d),2) 'deg']);
            imagesc(area.C_info_med(areaI,areaI,d) ./ area.est_variances_clean(areaI,d));
            axis image; colorbar;

            subplot(7,6,cnt+3); hold on;
            histogram(area.C_info_med(:,:,d) ./ area.est_variances_clean(:,d),linspace(-1,1,25));
            xlabel('correlation'); ylabel('log freq'); set(gca,'YScale','log');

            subplot(7,6,cnt+4); hold on; title(['strong corrs, disparity = ' num2str(cntr_disp(d),2) 'deg']);
            imagesc(area.C_info_strong(areaI,areaI,d) ./ area.est_variances_clean(areaI,d));
            axis image; colorbar;

            subplot(7,6,cnt+5); hold on;
            histogram(area.C_info_strong(:,:,d) ./ area.est_variances_clean(:,d),linspace(-1,1,25));
            xlabel('correlation'); ylabel('log freq'); set(gca,'YScale','log');

            cnt = cnt + 6;
        end
        saveas(figCorrs,['./plots/NoiseCorrelations/CorrelationMatrix_InfoLimiting_' this_area '.png']);
    end

    % remove the 3D matrices from structures because they inflate the file size
    area = rmfield(area,'C_info_weak');
    area = rmfield(area,'C_info_med');
    area = rmfield(area,'C_info_strong');

    % store results for each area in a separate matrix
    switch this_area
        case 'V1';      V1 = area;
        case 'V2';      V2 = area;
        case 'MT';      MT = area;
    end

end

%% plots

% confirm that information limiting correlations reduce overall FI

figure; hold on;
set(gcf,'Position',[10 10 700 200]);
title('information limiting correlations'); h = [];
h(1) = plot([1 2 3 4],[V1.FI_gauss_identity_total V1.FI_gauss_info_weak_total V1.FI_gauss_info_med_total V1.FI_gauss_info_strong_total],'o-','color',ColorIt('b'),'markerfacecolor',ColorIt('b'));
h(2) = plot([1 2 3 4],[V2.FI_gauss_identity_total V2.FI_gauss_info_weak_total V2.FI_gauss_info_med_total V2.FI_gauss_info_strong_total],'o-','color',ColorIt('g'),'markerfacecolor',ColorIt('g'));
h(3) = plot([1 2 3 4],[MT.FI_gauss_identity_total MT.FI_gauss_info_weak_total MT.FI_gauss_info_med_total MT.FI_gauss_info_strong_total],'o-','color',ColorIt('r'),'markerfacecolor',ColorIt('r'));
legend(h,'V1','V2','MT');
xticks([1 2 3 4]); box on; ylim([0 2e6]); xlim([0.75 4.25]); yticks([0 0.5*10^6 1*10^6 1.5*10^6 2*10^6])
xticklabels({'none','weak','medium','strong'}); ylabel('total Fisher Information')

saveas(gcf,['./plots/NoiseCorrelations/TotalInformation_limiting.png']);
saveas(gcf,['./plots/NoiseCorrelations/TotalInformation_limiting.svg']);


% plot the information limited population FIs

figure; hold on
plot(cntr_disp, V1.FI_gauss_info_weak,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_gauss_info_weak,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_gauss_info_weak,'color',ColorIt('r'),'linewidth',4);

title('Gaussian Noise, Weak Information Limiting Correlations')
axis square; box on;
xlabel('horizontal disparity (deg)');
ylabel('Fisher Information (normalized)');
set(gca,'fontsize',20); ylim([0 0.2]);
saveas(gcf,'./plots/NoiseCorrelations/FI_InfoWeakGaussian.png');
saveas(gcf,'./plots/NoiseCorrelations/FI_InfoWeakGaussian.svg');

figure; hold on
plot(cntr_disp, V1.FI_gauss_info_med,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_gauss_info_med,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_gauss_info_med,'color',ColorIt('r'),'linewidth',4);

title('Gaussian Noise, Medium Information Limiting Correlations')
axis square; box on;
xlabel('horizontal disparity (deg)');
ylabel('Fisher Information (normalized)');
set(gca,'fontsize',20); ylim([0 0.2]);
saveas(gcf,'./plots/NoiseCorrelations/FI_InfoMediumGaussian.png');
saveas(gcf,'./plots/NoiseCorrelations/FI_InfoMediumGaussian.svg');

figure; hold on
plot(cntr_disp, V1.FI_gauss_info_strong,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_gauss_info_strong,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_gauss_info_strong,'color',ColorIt('r'),'linewidth',4);

title('Gaussian Noise, Strong Information Limiting Correlations')
axis square; box on;
xlabel('horizontal disparity (deg)');
ylabel('Fisher Information (normalized)');
set(gca,'fontsize',20); ylim([0 0.2]);
saveas(gcf,'./plots/NoiseCorrelations/FI_InfoStrongGaussian.png');
saveas(gcf,'./plots/NoiseCorrelations/FI_InfoStrongGaussian.svg');

save('./analysisFiles/physio/results_populationFI_withCorrelations.mat',...
    'areas','V1','V2','MT','statsOmni','statsPairs','V1orig','V2orig','MTorig','subInds','cntr_disp','edges_disp');

