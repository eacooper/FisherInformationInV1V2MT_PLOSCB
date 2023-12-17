function [V1p,V2p,MTp] = permuteFits(V1,V2,MT,origin,target,parameter)

% Given a set of real fits to neurons from a cortical region, resample
% parameter fits given a f

% Define target and origin datasets
switch origin
    case 'V1'
        oriDat = V1;
    case 'V2'
        oriDat = V2;
    case 'MT'
        oriDat = MT;
end

switch target
    case 'V1'
        targDat = V1;
    case 'V2'
        targDat = V2;
    case 'MT'
        targDat = MT;
end

% Define parameter of interest
% Gabor params: [pedestal amplitude Gaussian_mean Gaussian_sigma frequency phase]
switch parameter
    case 'offset'
        pind = 1;
    case 'amplitude'
        pind = 2;
    case 'envMean'
        pind = 3;
    case 'envWidth'
        pind = 4;
    case 'frequency'
        pind = 5;
    case 'phase'
        pind = 6;
end

% Generate distribution for target area & parameter
theseParam    = targDat.P(:,pind);
numSamps      = size(oriDat.P,1);

[thisDens,xi] = ksdensity(theseParam);

% Resample from target distribution to create new origin distribution
pSamp = nan(numSamps,1);
for ii = 1:numSamps
    pSamp(ii) = pinky(xi,1,thisDens);
end

% Save to new permuted dataset
V1p = V1;
V2p = V2;
MTp = MT;

switch origin
    case 'V1'
        V1p.P(:,pind) = pSamp;
    case 'V2'
        V2p.P(:,pind) = pSamp;
    case 'MT'
        MT.P(:,pind)  = pSamp;
end

end