function area = fitNeuronalVariances(area_title,area,plotallfits)

minReps = 3;

for n = 1:size(area.allresp,1)

    % grab mean, variance, and repeats
    these_means = area.experiments{n}.dat(:,2);
    these_vars = area.experiments{n}.dat(:,3);
    these_reps = area.experiments{n}.dat(:,4);

    % only use disparities with at least 3 repeats
    these_means = these_means(these_reps >= minReps);
    these_vars = these_vars(these_reps >= minReps);

    % fit a line to means to predict variances
    X = [ones(size(these_means)) these_means];
    %b = regress(these_vars,X);

    % fit a line to means to predict variances, with the constraint that
    % the variance must be greater than zero
    options = optimset('Algorithm','interior-point','Display', 'off');
    B=lsqlin(X,these_vars,[],[],[],[],[0,-inf],[],[],options);

    area.mean_var_intercept(n) = B(1);
    area.mean_var_slope(n)     = B(2);

    % store for plotting/analysis
    all_means{n} = these_means;
    all_vars{n} = these_vars;

end


if(plotallfits)

    thisfig = figure; hold on;

    % counters
    pcnt = 1;
    fcnt = 1;
    for n = 1:length(area.experiments)

        if mod(n,66) == 1
            if fcnt > 1
                saveas(gcf,['./plots/AssessFits/AllFits/MeanVariance_' area_title '_' num2str(fcnt-1) '.png']);
                close gcf;
            end
            thisfig = figure; hold on;
            setupfig(18,10,10);
            pcnt = 1;
            fcnt = fcnt + 1;
        end

        subplot(6,11,pcnt); hold on; title([ num2str(n) ]);

        % data
        scatter(all_means{n},all_vars{n},'k','filled');

        % our fit
        plot([min(all_means{n}) max(all_means{n})],area.mean_var_slope(n)*[min(all_means{n}) max(all_means{n})] + area.mean_var_intercept(n),'b-');

        % Poisson identity line for reference
        plot([min(all_means{n}) max(all_means{n})],[min(all_means{n}) max(all_means{n})],'r--');

        box on;

        xlabel('spike rate mean (sps)'); ylabel('spike rate variance (sps2)');
        pcnt = pcnt + 1;

    end
    saveas(gcf,['./plots/AssessFits/AllFits/MeanVariance_' area_title '_' num2str(fcnt-1) '.png']);

    if exist('thisfig')
        close(thisfig);
    end
end