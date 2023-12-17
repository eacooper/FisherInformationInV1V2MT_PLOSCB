% This script finds the best power law that relates the measured disparity
% distributions from the BORIS image dataset to the FI distributions from
% the neural datasets

clear all; close all;

% Load in FI data
load('./analysisFiles/physio/results_populationFI_withCorrelations.mat');

% Tasks for disparity statistics
tasks = {'Walk','Sando'};

% we'll make a plot of the Error for each power to make sure each one
% reaches a minimum within the range we're testing
figError = figure;

% Loop over each area
for a = 1:length(areas)

    % Load in disparity histograms from BORIS image set, sampled to match this area
    sndo = load(['./analysisFiles/disparityStats/dispHist' areas{a} '_sando.mat']);
    wlk  = load(['./analysisFiles/disparityStats/dispHist' areas{a} '_walking.mat']);

    % grab the physio data and disparity histogram associated with each area
    switch areas{a}

        case 'V1'
            area = V1;
            Sando(a,:)      = sndo.dispHistV1./sum(sndo.dispHistV1);    % sando scene stats
            Walk(a,:)       = wlk.dispHistV1./sum(wlk.dispHistV1);      % walking scene stats
            clr = 'b';
        case 'V2'
            area = V2;
            Sando(a,:)  = sndo.dispHistV2./sum(sndo.dispHistV2);
            Walk(a,:)   = wlk.dispHistV2./sum(wlk.dispHistV2);
            clr = 'g';
        case 'MT'
            area = MT;
            Sando(a,:)  = sndo.dispHistMT./sum(sndo.dispHistMT);
            Walk(a,:)   = wlk.dispHistMT./sum(wlk.dispHistMT);
            clr = 'r';
        otherwise;  error('invalid area')

    end


    % Calculate best fitting power law using the neurally-defined KSDs

    % Define powers to test
    ps = linspace(0.1,4,300);

    % Loop over set of powers
    for p = 1:length(ps)

        Sando_pow = (Sando(a,:).^ps(p))./sum(Sando(a,:).^ps(p));
        Walk_pow = (Walk(a,:).^ps(p))./sum(Walk(a,:).^ps(p));

        % error for main Poisson statistics analysis
        Sando_err(p) = mean(abs(Sando_pow-area.FI_poiss_pop));
        Walk_err(p) = mean(abs(Walk_pow-area.FI_poiss_pop));

        % gaussian noise with independence
        Sando_err_gauss_var(p) = mean(abs(Sando_pow-area.FI_gauss_var_pop));
        Walk_err_gauss_var(p) = mean(abs(Walk_pow-area.FI_gauss_var_pop));

        % ... or with information-limiting correlations
        Sando_err_gauss_info_weak(p) = mean(abs(Sando_pow-area.FI_gauss_info_weak));
        Walk_err_gauss_info_weak(p) = mean(abs(Walk_pow-area.FI_gauss_info_weak));

        Sando_err_gauss_info_med(p) = mean(abs(Sando_pow-area.FI_gauss_info_med));
        Walk_err_gauss_info_med(p) = mean(abs(Walk_pow-area.FI_gauss_info_med));

        Sando_err_gauss_info_strong(p) = mean(abs(Sando_pow-area.FI_gauss_info_strong));
        Walk_err_gauss_info_strong(p) = mean(abs(Walk_pow-area.FI_gauss_info_strong));


    end

    % Find power that most closely matches FI distribution and disparity stats

    % Poisson
    Sando_err_min(a) = ps((Sando_err == min(Sando_err)));
    Walk_err_min(a)  = ps((Walk_err == min(Walk_err)));

    % Gaussians
    Sando_err_min_gauss(a,1) = ps((Sando_err_gauss_var == min(Sando_err_gauss_var)));
    Walk_err_min_gauss(a,1)  = ps((Walk_err_gauss_var == min(Walk_err_gauss_var)));
    Sando_err_min_gauss(a,2) = ps((Sando_err_gauss_info_weak == min(Sando_err_gauss_info_weak)));
    Walk_err_min_gauss(a,2)  = ps((Walk_err_gauss_info_weak == min(Walk_err_gauss_info_weak)));
    Sando_err_min_gauss(a,3) = ps((Sando_err_gauss_info_med == min(Sando_err_gauss_info_med)));
    Walk_err_min_gauss(a,3)  = ps((Walk_err_gauss_info_med == min(Walk_err_gauss_info_med)));
    Sando_err_min_gauss(a,4) = ps((Sando_err_gauss_info_strong == min(Sando_err_gauss_info_strong)));
    Walk_err_min_gauss(a,4)  = ps((Walk_err_gauss_info_strong == min(Walk_err_gauss_info_strong)));

    % Plot error function
    figure(figError);
    hold on;

    plot(ps,Sando_err,'-','color',ColorIt(clr),'linewidth',2);
    plot(ps,Walk_err,'--','color',ColorIt(clr),'linewidth',2);

    % Repeat Poisson analysis with each of the bootstrapped disparity samples
    lfi = area.FI_poiss_boot;

    % number bootstraps
    nboots = size(lfi,1);

    % for each boot
    for n = 1:nboots

        % for each power
        for p = 1:length(ps)

            Sando_pow = (Sando(a,:).^ps(p))./sum(Sando(a,:).^ps(p));
            Sando_err(p) = mean(abs(Sando_pow-lfi(n,:)));

            Walk_pow = (Walk(a,:).^ps(p))./sum(Walk(a,:).^ps(p));
            Walk_err(p) = mean(abs(Walk_pow-lfi(n,:)));

        end

        % find minimum
        Sando_err_min_boot(a,n) = ps((Sando_err == min(Sando_err)));
        Walk_err_min_boot(a,n) = ps((Walk_err == min(Walk_err)));

    end

    switch areas{a}

        case 'V1';  V1 = area;
        case 'V2';  V2 = area;
        case 'MT';  MT = area;
        otherwise;  error('invalid area')

    end

end


% Add legend and save Error plot
figure(figError);
hold on;
legend('Sando v1','Walk v1','Sando v2','Walk v2','Sando mt','Walk mt','location','northeastoutside');
title('Error with FI_poiss');
xlabel('power'); ylabel('Error')
saveas(figError,'./plots/NDSComparison/example_gridsearch_error.svg');

%---------------------------------------------------------%
%% Plot a histogram of the best fit powers for Poisson, fit a Gaussian to error distribution,
% then determine power CIs based on that Gaussian

% histogram bin placement
bins = linspace(0.25,2,35);

for t = 1:2

    task = tasks{t};

    switch task
        case 'Walk'
            err_min = Walk_err_min;
            err_min_boot = Walk_err_min_boot;

        case 'Sando'
            err_min = Sando_err_min;
            err_min_boot = Sando_err_min_boot;

    end

    % Plot histograms of best fit powers
    f1{t} = figure;
    f1{t}.Position = [500+500*(t-1) 100 650 400];
    hold on;

    h1 = histogram(err_min_boot(1,:),bins,'Normalization','pdf','FaceAlpha',0.5);
    h2 = histogram(err_min_boot(2,:),bins,'Normalization','pdf','FaceAlpha',0.5);
    h3 = histogram(err_min_boot(3,:),bins,'Normalization','pdf','FaceAlpha',0.5);

    h1.FaceColor = ColorIt('b');
    h2.FaceColor = ColorIt('g');
    h3.FaceColor = ColorIt('r');

    box on;

    % Fit Gaussians to power histograms and add to plot
    plt_points = linspace(0.2,2,200);

    [muHat_V1,sigmaHat_V1,muCI_V1,sigmaCI_V1] = normfit(err_min_boot(1,:));
    plot(plt_points,normpdf(plt_points,muHat_V1,sigmaHat_V1),'color',ColorIt('b'),'linewidth',3);

    [muHat_V2,sigmaHat_V2,muCI_V2,sigmaCI_V2] = normfit(err_min_boot(2,:));
    plot(plt_points,normpdf(plt_points,muHat_V2,sigmaHat_V2),'color',ColorIt('g'),'linewidth',3);

    [muHat_MT,sigmaHat_MT,muCI_MT,sigmaCI_MT] = normfit(err_min_boot(3,:));
    plot(plt_points,normpdf(plt_points,muHat_MT,sigmaHat_MT),'color',ColorIt('r'),'linewidth',3);

    xlim([min(bins) max(bins)]);
    set(gca,'fontsize',20,'xlim',[0.2 2],'ylim',[0 20]);
    xlabel('Power');
    ylabel('Bootstrap count');
    title(task);

    % Cohen's D
    D_V1_MT = (muHat_V1 - muHat_MT) / sqrt( ( (nboots-1)*(sigmaHat_V1^2) + (nboots-1)*(sigmaHat_MT^2) ) / (2*nboots - 2) );
    D_V2_MT = (muHat_V2 - muHat_MT) / sqrt( ( (nboots-1)*(sigmaHat_V2^2) + (nboots-1)*(sigmaHat_MT^2) ) / (2*nboots - 2) );
    D_V1_V2 = (muHat_V1 - muHat_V2) / sqrt( ( (nboots-1)*(sigmaHat_V1^2) + (nboots-1)*(sigmaHat_V2^2) ) / (2*nboots - 2) );

    % report V1/MT and V2/MT ratios of the best fit power laws for Poisson
    ratio_V1_MT = err_min(1)/err_min(3);
    ratio_V2_MT = err_min(2)/err_min(3);
    ratio_V1_V2 = err_min(1)/err_min(2);

    % put mus, Ds and ratios together
    stats = [err_min(1) err_min(2) err_min(3) D_V1_MT D_V2_MT D_V1_V2 ratio_V1_MT ratio_V2_MT ratio_V1_V2];

    switch task
        case 'Walk'
            statsWalk = stats;

        case 'Sando'
            statsSando = stats;

    end

end

% put stats into tables to save and print to workspace
% Define the values for each row
analysis = {'best V1', 'best V2', 'best MT', 'D V1-MT', 'D V2-MT','D V1-V2', 'ratio V1-MT','ratio V2-MT', 'ratio V1-V2'};

% Create the table
statsTasks = table(analysis', statsWalk', statsSando','VariableNames', {'analysis', 'walk', 'sando'});

% Display the table
disp('Mus, Cohens Ds, and ratios')
disp(statsTasks);

%---------------------------------------------------------%
% Overlay FI with best fit power

% Food preparation
f2 = figure;
f2.Position = [100 700 575 570];
hold on;

h(1) = plot(cntr_disp, V1.FI_poiss_pop,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_poiss_pop,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_poiss_pop,'color',ColorIt('r'),'linewidth',4);

h(2) = plot(cntr_disp, (Sando(1,:).^Sando_err_min(1))./sum(Sando(1,:).^Sando_err_min(1)),':','color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, (Sando(2,:).^Sando_err_min(2))./sum(Sando(2,:).^Sando_err_min(2)),':','color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, (Sando(3,:).^Sando_err_min(3))./sum(Sando(3,:).^Sando_err_min(3)),':','color',ColorIt('r'),'linewidth',4);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)');
ylabel('Probability density');
set(gca,'fontsize',20);
legend(h,'Pop. FI','Disp. Dist.^p');
title('sando');


% Navigation
f3 = figure;
f3.Position = [100 100 575 570];
hold on;

h(1) = plot(cntr_disp, V1.FI_poiss_pop,'color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, V2.FI_poiss_pop,'color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, MT.FI_poiss_pop,'color',ColorIt('r'),'linewidth',4);

h(2) = plot(cntr_disp, (Walk(1,:).^Walk_err_min(1))./sum(Walk(1,:).^Walk_err_min(1)),':','color',ColorIt('b'),'linewidth',4);
plot(cntr_disp, (Walk(2,:).^Walk_err_min(2))./sum(Walk(2,:).^Walk_err_min(2)),':','color',ColorIt('g'),'linewidth',4);
plot(cntr_disp, (Walk(3,:).^Walk_err_min(3))./sum(Walk(3,:).^Walk_err_min(3)),':','color',ColorIt('r'),'linewidth',4);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)');
ylabel('Probability density');
set(gca,'fontsize',20);
legend(h,'Pop. FI','Disp. Dist.^p');
title('walk');


% plot V1/MT and V2/MT ratios of the best fit power laws for each noise model
Walk_ratio_V1_MT = statsWalk(7);
Walk_ratio_V2_MT = statsWalk(8);
Walk_ratio_V1_V2 = statsWalk(9);

Sando_ratio_V1_MT = statsSando(7);
Sando_ratio_V2_MT = statsSando(8);
Sando_ratio_V1_V2 = statsSando(9);

Walk_ratio_V1_MT_gauss = Walk_err_min_gauss(1,:)./ Walk_err_min_gauss(3,:);
Walk_ratio_V2_MT_gauss = Walk_err_min_gauss(2,:)./ Walk_err_min_gauss(3,:);
Walk_ratio_V1_V2_gauss = Walk_err_min_gauss(1,:)./ Walk_err_min_gauss(2,:);

Sando_ratio_V1_MT_gauss = Sando_err_min_gauss(1,:)./ Sando_err_min_gauss(3,:);
Sando_ratio_V2_MT_gauss = Sando_err_min_gauss(2,:)./ Sando_err_min_gauss(3,:);
Sando_ratio_V1_V2_gauss = Sando_err_min_gauss(1,:)./ Sando_err_min_gauss(2,:);

f4 = figure;
f4.Position = [100   454   746   216];
hold on;

xvals = [1.5 2 2.5];
for m = 1:4

        hR(1) = plot(xvals,[Walk_ratio_V1_MT_gauss(m) Walk_ratio_V2_MT_gauss(m) Walk_ratio_V1_V2_gauss(m)],'o-','color',ColorIt('o'),'markerfacecolor',ColorIt('o'));
        hR(2) = plot(xvals,[Sando_ratio_V1_MT_gauss(m) Sando_ratio_V2_MT_gauss(m) Sando_ratio_V1_V2_gauss(m)],'s-','color',ColorIt('m'),'markerfacecolor',ColorIt('m'));

        % plot Poisson for comparison
        plot(xvals,[Walk_ratio_V1_MT Walk_ratio_V2_MT Walk_ratio_V1_V2],':','color',ColorIt('o'));
        plot(xvals,[Sando_ratio_V1_MT Sando_ratio_V2_MT Sando_ratio_V1_V2],':','color',ColorIt('m'));

        if m == 1; text(xvals(2), 1.8, 'Ind. Gaussian (fitted)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif m == 2; text(xvals(2), 1.9, 'Weak Corr.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif m == 3; text(xvals(2), 1.9, 'Medium Corr.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        elseif m == 4; text(xvals(2), 1.9, 'Strong Corr.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end

    xvals = xvals + 3;

end

plot([0 13],[1 1],'k-');
xticks([1.5 2 2.5 4.5 5 5.5 7.5 8 8.5 10.5 11 11.5]); xticklabels({'V1/MT','V2/MT','V1/V2','V1/MT','V2/MT','V1/V2','V1/MT','V2/MT','V1/V2','V1/MT','V2/MT','V1/V2'});
xtickangle(45);
ylabel('power law ratio');
ylim([0.5 2]); xlim([0.5 12.5]); yticks([0.5 1 1.5 2]);
legend(hR,'navigation','food preparation');

% f1: Histograms and best fit Gaussians to histogram
saveas(f1{1},'./plots/NDSComparison/power_Walk_histogram.svg');
saveas(f1{2},'./plots/NDSComparison/power_Sando_histogram.svg');

% f2: Best fit power law xformed disparity distribution vs. FI (food prep)
saveas(f2,'./plots/NDSComparison/power_Sando_fits.svg');

% f3: Best fit power law xformed disparity distribution vs. FI (navigation)
saveas(f3,'./plots/NDSComparison/power_Walk_fits.svg');

saveas(f4,'./plots/NDSComparison/power_Ratios.svg');

save('./analysisFiles/results_populationFI_withCorrelations_NDS_comparison.mat',...
    'areas','V1','V2','MT','statsOmni','statsPairs','statsTasks','V1orig','V2orig','MTorig','subInds','cntr_disp','edges_disp');



