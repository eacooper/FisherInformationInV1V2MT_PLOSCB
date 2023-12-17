function [h] = shadeplot(ctIn,ebIn,hIn,color_in,alpha,linewidth)
    % Plots a filled shaded region instead of error bars or IQRs. Assumes
    % input arguments are at most nx2 matrices, lower CI/IQR is the first dim.)
    %
    % Usage: shadeplot(centTendencyMet,ErrorBarMet,hAxis,colorIn,alpha,lineWidth)

    % Make sure input data are all oriented with longest dim first
    ctSz = size(ctIn);
    ebSz = size(ebIn);
    hSz  = size(hIn);

    if ctSz(2) > ctSz(1)
        ctIn = ctIn';
    end

    if ebSz(2) > ebSz(1)
        ebIn = ebIn';
    end

    if hSz(2) > hSz(1)
        hIn = hIn';
    end

    % Check to see if mean/standard deviation do not change with horizontal 
    if isscalar(ctIn)
        F = ctIn.*ones(horizIn,1);
    else
        F = ctIn;
    end
    
    if isscalar(ebIn)
        F_std = ebIn.*ones(horizIn,1);
    else
        F_std = ebIn;
    end

    % Plot std or IQR as a filled polygon
    if size(ebIn,2)>1
        % If there are asymmetric errorbars
        xVals = [hIn' fliplr(hIn') hIn(1)];
        yVals = [ebIn(:,1)' fliplr(ebIn(:,2)') ebIn(1,1)];
    else
        % If errorbars are symmetric
        xVals = [hIn' fliplr(hIn')];
        yVals = [F + F_std fliplr(F - F_std)]';
    end

    fill(xVals,yVals,color_in,'FaceAlpha',alpha,'linestyle','none');

    % Plot central tendency metric
    h = plot(hIn,ctIn,'Color',color_in,'lineWidth',linewidth);

end