function getDispStats(imSet,horizontal_disparity,numSamps,ecc,edges_disp,res,resampIter,V1densityMat,V2densityMat,MTdensityMat,circDensityMat)
% For image set what are the disparity statistics?

% Collect disparities within masks and histogram them

% initialize
numIms = size(horizontal_disparity,3);

dispHistV1 = nan(numIms,res-1);
dispHistV2 = nan(numIms,res-1);
dispHistMT = nan(numIms,res-1);
dispHistCirc = nan(numIms,res-1);

disp(['Iteration ' num2str(resampIter) ' for ' imSet])
% Loop over all images in dataset
for ii = 1:numIms

    if mod(ii,1000) == 0
        disp(['Running image ',num2str(ii),'/',num2str(numIms)]);
    end

    % For bootstrapping, select a random image from set
    if resampIter ~= 0
        imInd = randi(numIms);
    else
        imInd = ii;
    end

    % Make a mask of the visual field subregion of interest
    mask = ecc<10;

    % BORIS image set has row 1 = -10.3deg, row 2 = -10.3deg, etc so let's
    % flip it to match RF probability densities
    dispCrop = flipud(squeeze(horizontal_disparity(:,:,imInd)));

    % Resample disparities based on RF location probability densities
    imSize   = size(MTdensityMat,1);

    xSampV1    = nan(numSamps,1);
    ySampV1    = nan(numSamps,1);
    xSampV2    = nan(numSamps,1);
    ySampV2    = nan(numSamps,1);
    xSampMT    = nan(numSamps,1);
    ySampMT    = nan(numSamps,1);
    xSampCirc  = nan(numSamps,1);
    ySampCirc  = nan(numSamps,1);

    % First modify KSD mats so we mask out pixels with undefined
    % disparities and/or unwanted ROIs in the VF
    undefImMask = ~isnan(dispCrop);

    V1densityMatMasked   = V1densityMat.*mask.*undefImMask;
    V2densityMatMasked   = V2densityMat.*mask.*undefImMask;
    MTdensityMatMasked   = MTdensityMat.*mask.*undefImMask;
    circDensityMatMasked = circDensityMat.*mask.*undefImMask;

    % Sometimes the KSD plots are all zero after masking if some of the
    % images are full of undefined regions. If one of these images is
    % encountered, just skip it
    check(1) = sum(V1densityMatMasked(:));
    check(2) = sum(V2densityMatMasked(:));
    check(3) = sum(MTdensityMatMasked(:));
    check(4) = sum(circDensityMatMasked(:));

    allCheck = sum([check(1) == 0; check(2) == 0; check(3) == 0; check(4) == 0]);

    if allCheck
        continue
    end

    % Sample indices of disparity plots based on KSD
    for jj = 1:numSamps

        [ySampV1(jj),xSampV1(jj)]     = pinky(1:imSize,1:imSize,V1densityMatMasked);
        [ySampV2(jj),xSampV2(jj)]     = pinky(1:imSize,1:imSize,V2densityMatMasked);
        [ySampMT(jj),xSampMT(jj)]     = pinky(1:imSize,1:imSize,MTdensityMatMasked);
        [ySampCirc(jj),xSampCirc(jj)] = pinky(1:imSize,1:imSize,circDensityMatMasked);

    end

    % convert to linear indices
    sampIndV1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1,ySampV1);
    sampIndV2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV2,ySampV2);
    sampIndMT   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampMT,ySampMT);
    sampIndCirc = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampCirc,ySampCirc);

    % Select disparities using these indices
    dispResampV1   = dispCrop(sampIndV1);
    dispResampV2   = dispCrop(sampIndV2);
    dispResampMT   = dispCrop(sampIndMT);
    dispResampCirc = dispCrop(sampIndCirc);

    % Calculate final histograms
    dispHistV1(ii,:)   = histcounts(dispResampV1,edges_disp);
    dispHistV2(ii,:)   = histcounts(dispResampV2,edges_disp);
    dispHistMT(ii,:)   = histcounts(dispResampMT,edges_disp);
    dispHistCirc(ii,:) = histcounts(dispResampCirc,edges_disp);
end


% Collapse across images
dispHistV1   = sum(dispHistV1,1,'omitnan');
dispHistV2   = sum(dispHistV2,1,'omitnan');
dispHistMT   = sum(dispHistMT,1,'omitnan');
dispHistCirc = sum(dispHistCirc,1,'omitnan');

% Normalize to get PDF
dispHistV1   = dispHistV1/(sum(dispHistV1)*diff(edges_disp(1:2)));
dispHistV2   = dispHistV2/(sum(dispHistV2)*diff(edges_disp(1:2)));
dispHistMT   = dispHistMT/(sum(dispHistMT)*diff(edges_disp(1:2)));
dispHistCirc = dispHistCirc/(sum(dispHistCirc)*diff(edges_disp(1:2)));


%% Save
if resampIter ~= 0
    % bootstrap
    suffix = num2str(resampIter);

    save(['./analysisFiles/disparityStats/bootstraps/dispHistV1_',imSet,suffix],'dispHistV1');
    save(['./analysisFiles/disparityStats/bootstraps/dispHistV2_',imSet,suffix],'dispHistV2');
    save(['./analysisFiles/disparityStats/bootstraps/dispHistMT_',imSet,suffix],'dispHistMT');
    save(['./analysisFiles/disparityStats/bootstraps/dispHistCirc_',imSet,suffix],'dispHistCirc');
else
    % main data
    suffix = '';

    save(['./analysisFiles/disparityStats/dispHistV1_',imSet,suffix],'dispHistV1');
    save(['./analysisFiles/disparityStats/dispHistV2_',imSet,suffix],'dispHistV2');
    save(['./analysisFiles/disparityStats/dispHistMT_',imSet,suffix],'dispHistMT');
    save(['./analysisFiles/disparityStats/dispHistCirc_',imSet,suffix],'dispHistCirc');
end




end


