%% In the main analysis, we computed fisher information for each neuron
% assuming Poisson statistics, then summed them up because we assumed
% statistical independence. In this control analysis, we test how the
% results might change if the neurons has statistical dependences. This is
% a supplemental script that explores local correlations, rather than
% information limiting correlations

clear all; close all;

% Load Gabor fitting results of the neurons that we want to analyze
load('analysisFiles/physio/results_populationFI.mat');

% load symbolic equations for Gabor and first derivate
define_Gabor_tuning_curve;

% Now, we define a few pairwise correlation matrices for local correlations

% specify the amplitude and sigma of a Gaussian to tell us how correlated two neurons are as a
% function of how different their disparity preference is
max_weak = 0.1;
sig_weak = 0.5;

max_med = 0.2;
sig_med = 0.5;

max_strong = 0.3;
sig_strong = 0.5;

% for each area
for a = 1:length(areas)

    this_area = areas{a};

    % grab the data structure for this area
    switch this_area
        case 'V1';      area = V1;
        case 'V2';      area = V2;
        case 'MT';      area = MT;
    end

    display(['making correlation matrices for ' this_area]);

    % how many neurons in this area?
    numNeurons = size(area.allresp,1);

    % for each pair of neurons
    for i = 1:numNeurons
        for j = 1:numNeurons

            % difference in preferred disparity
            disp_pref_diff = area.pref_disp(i) - area.pref_disp(j);

            % assert desired correlations
            if i == j
                area.C_identity_corr(i,j) = 1;
                area.C_local_weak_corr(i,j) = 1;
                area.C_local_med_corr(i,j) = 1;
                area.C_local_strong_corr(i,j) = 1;
            else
                area.C_identity_corr(i,j) = 0;
                area.C_local_weak_corr(i,j) = max_weak*exp(-(disp_pref_diff^2)/(sig_weak^2));
                area.C_local_med_corr(i,j) = max_med*exp(-(disp_pref_diff^2)/(sig_med^2));
                area.C_local_strong_corr(i,j) = max_strong*exp(-(disp_pref_diff^2)/(sig_strong^2));
            end


            % for each disparity, convert correlations to covariance by scaling by the std's of these two neurons
            for d = 1:length(cntr_disp)

                area.C_identity(i,j,d) = area.C_identity_corr(i,j)*sqrt(area.est_variances_clean(i,d))*sqrt(area.est_variances_clean(j,d));
                area.C_local_weak(i,j,d) = area.C_local_weak_corr(i,j)*sqrt(area.est_variances_clean(i,d))*sqrt(area.est_variances_clean(j,d));
                area.C_local_med(i,j,d) = area.C_local_med_corr(i,j)*sqrt(area.est_variances_clean(i,d))*sqrt(area.est_variances_clean(j,d));
                area.C_local_strong(i,j,d) = area.C_local_strong_corr(i,j)*sqrt(area.est_variances_clean(i,d))*sqrt(area.est_variances_clean(j,d));

            end

        end

    end

    switch this_area
        case 'V1';      V1 = area;
        case 'V2';      V2 = area;
        case 'MT';      MT = area;
    end

end

%% plots

% plot functions that determine weak, medium and strong local correlations
figure; hold on; h = [];
h(1) = plot(cntr_disp,max_weak*exp(-(cntr_disp.^2)./(sig_weak^2)),'k:');
h(2) = plot(cntr_disp,max_med*exp(-(cntr_disp.^2)./(sig_med^2)),'k--');
h(3) = plot(cntr_disp,max_strong*exp(-(cntr_disp.^2)./(sig_strong^2)),'k-');
legend(h,'weak correlations','medium correlations','strong correlations');
xlabel('difference in disparity pref (deg)'); ylabel('correlation');
saveas(gcf,'./plots/NoiseCorrelations/Local/CorrelationMatrix_LocalCorrFunctions.png');

% sort neurons by their preferred disparity
[~,V1I] = sort(V1.pref_disp);
[~,V2I] = sort(V2.pref_disp);
[~,MTI] = sort(MT.pref_disp);

% weak local correlations
figure; hold on; set(gcf,'Position',[10 10 1000 1000]); sgtitle('Weak local correlations')
subplot(1,3,1); hold on; title('V1'); imagesc(V1.C_local_weak_corr(V1I,V1I)); axis image; colorbar;
subplot(1,3,2); hold on; title('V2'); imagesc(V2.C_local_weak_corr(V2I,V2I)); axis image; colorbar;
subplot(1,3,3); hold on; title('MT'); imagesc(MT.C_local_weak_corr(MTI,MTI)); axis image; colorbar;
saveas(gcf,'./plots/NoiseCorrelations/Local/CorrelationMatrix_LocalWeak.png');

% medium
figure; hold on; set(gcf,'Position',[10 10 1000 1000]); sgtitle('Medium local correlations')
subplot(1,3,1); hold on; title('V1'); imagesc(V1.C_local_med_corr(V1I,V1I)); axis image; colorbar;
subplot(1,3,2); hold on; title('V2'); imagesc(V2.C_local_med_corr(V2I,V2I)); axis image; colorbar;
subplot(1,3,3); hold on; title('MT'); imagesc(MT.C_local_med_corr(MTI,MTI)); axis image; colorbar;
saveas(gcf,'./plots/NoiseCorrelations/Local/CorrelationMatrix_LocalMedium.png');

% strong
figure; hold on; set(gcf,'Position',[10 10 1000 1000]); sgtitle('Strong local correlations')
subplot(1,3,1); hold on; title('V1'); imagesc(V1.C_local_strong_corr(V1I,V1I)); axis image; colorbar;
subplot(1,3,2); hold on; title('V2'); imagesc(V2.C_local_strong_corr(V2I,V2I)); axis image; colorbar;
subplot(1,3,3); hold on; title('MT'); imagesc(MT.C_local_strong_corr(MTI,MTI)); axis image; colorbar;
saveas(gcf,'./plots/NoiseCorrelations/Local/CorrelationMatrix_LocalStrong.png');


%% How does the Gaussian FI change if we apply these correlation matrices?
% for each disparity, we take the product of a vector containing the first
% derivative of each neuron at that disparity transposted, times the correlation
% matrix, times the derivative vector untransposed

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

    % weak local correlations
    for d = 1:length(cntr_disp)
        thisC = area.C_local_weak(:,:,d);
        invC = inv(thisC);
        area.FI_gauss_local_weak(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);
    end
    area.FI_gauss_local_weak_total = sum(area.FI_gauss_local_weak);
    area.FI_gauss_local_weak =  area.FI_gauss_local_weak / area.FI_gauss_local_weak_total;

    % medium local correlations
    for d = 1:length(cntr_disp)
        thisC = area.C_local_med(:,:,d);
        invC = inv(thisC);
        area.FI_gauss_local_med(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);
    end
    area.FI_gauss_local_med_total = sum(area.FI_gauss_local_med);
    area.FI_gauss_local_med =  area.FI_gauss_local_med / area.FI_gauss_local_med_total;

    % strong local correlations
    for d = 1:length(cntr_disp)
        thisC = area.C_local_strong(:,:,d);
        invC = inv(thisC);
        area.FI_gauss_local_strong(d) = area.first_derivs(:,d)'*invC*area.first_derivs(:,d);
    end
    area.FI_gauss_local_strong_total = sum(area.FI_gauss_local_strong);
    area.FI_gauss_local_strong =  area.FI_gauss_local_strong / area.FI_gauss_local_strong_total;

    % store results for each area in a separate matrix
    switch this_area
        case 'V1';      V1 = area;
        case 'V2';      V2 = area;
        case 'MT';      MT = area;
    end

end

%% plot the various population FIs
figure; hold on
plot(cntr_disp, V1.FI_gauss_local_weak,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_gauss_local_weak,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_gauss_local_weak,'color',ColorIt('r'),'linewidth',4);

title('Gaussian Noise, Weak Local Correlations')
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Fisher Information (normalized)');
set(gca,'fontsize',20);
saveas(gcf,'./plots/NoiseCorrelations/Local/FI_LocalWeakGaussian.png');

figure; hold on
plot(cntr_disp, V1.FI_gauss_local_med,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_gauss_local_med,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_gauss_local_med,'color',ColorIt('r'),'linewidth',4);

title('Gaussian Noise, Medium Local Correlations')
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Fisher Information (normalized)');
set(gca,'fontsize',20);
saveas(gcf,'./plots/NoiseCorrelations/Local/FI_LocalMediumGaussian.png');

figure; hold on
plot(cntr_disp, V1.FI_gauss_local_strong,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_gauss_local_strong,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_gauss_local_strong,'color',ColorIt('r'),'linewidth',4);

title('Gaussian Noise, Strong Local Correlations')
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Fisher Information (normalized)');
set(gca,'fontsize',20);
saveas(gcf,'./plots/NoiseCorrelations/Local/FI_LocalStrongGaussian.png');


% plots

% confirm that local correlations should not affect overall FI much
figure; hold on;
set(gcf,'Position',[10 10 700 200]);
title('local correlations'); h = [];
h(1) = plot([1 2 3 4],[V1.FI_gauss_identity_total V1.FI_gauss_local_weak_total V1.FI_gauss_local_med_total V1.FI_gauss_local_strong_total],'o-','color',ColorIt('b'),'markerfacecolor',ColorIt('b'));
h(2) = plot([1 2 3 4],[V2.FI_gauss_identity_total V2.FI_gauss_local_weak_total V2.FI_gauss_local_med_total V2.FI_gauss_local_strong_total],'o-','color',ColorIt('g'),'markerfacecolor',ColorIt('g'));
h(3) = plot([1 2 3 4],[MT.FI_gauss_identity_total MT.FI_gauss_local_weak_total MT.FI_gauss_local_med_total MT.FI_gauss_local_strong_total],'o-','color',ColorIt('r'),'markerfacecolor',ColorIt('r'));
legend(h,'V1','V2','MT');
xticks([1 2 3 4]); box on; ylim([0 2e6]); xlim([0.75 4.25]); yticks([0 0.5*10^6 1*10^6 1.5*10^6 2*10^6])
xticklabels({'none','weak','medium','strong'}); ylabel('total Fisher Information')

saveas(gcf,['./plots/NoiseCorrelations/Local/TotalInformation_local.png']);
saveas(gcf,['./plots/NoiseCorrelations/Local/TotalInformation_local.svg']);
