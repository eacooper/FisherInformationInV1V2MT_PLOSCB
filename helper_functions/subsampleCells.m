function [subInds] = subsampleCells(V1,V2,MT)

% Take full sample of V1/V2/MT cells and resample them so we're only
% looking at those well-described by Gabor tuning, those with RF centers
% within 10deg of fixation, and roughly all in the same visual field
% regions

%% Only keep neurons with R2 > 0.75
r2Inds.V1 = V1.r2 > 0.75;
r2Inds.V2 = V2.r2 > 0.75;
r2Inds.MT = MT.r2 > 0.75;

display(['R2 .75 V1 =  ' num2str(sum(r2Inds.V1))]);
display(['R2 .75 V2 =  ' num2str(sum(r2Inds.V2))]);
display(['R2 .75 MT =  ' num2str(sum(r2Inds.MT))]);

%% remove cells from V1 and V2 where the fixation point was not the center
% of the screen, because in this case the conversion to disparity is
% incorrect
for x = 1:size(V1.experiments,2)
    fpV1 = V1.experiments{x}.alldata.fp;
    if all(isnan(fpV1)) || all(fpV1 == 0)
        fpInds.V1(x) = 1;
    else
        fpInds.V1(x) = 0;
    end
end

for x = 1:size(V2.experiments,2)
    fpV2 = V2.experiments{x}.alldata.fp;
    if all(isnan(fpV2)) || all(fpV2 == 0)
        fpInds.V2(x) = 1;
    else
        fpInds.V2(x) = 0;
    end
end

subInds.V1 = r2Inds.V1 & fpInds.V1;
subInds.V2 = r2Inds.V2 & fpInds.V2;
subInds.MT = r2Inds.MT;

display(['V1 fp 0:  ' num2str(sum(subInds.V1))]);
display(['V2 fp 0:  ' num2str(sum(subInds.V2))]);
display(['MT fp 0:  ' num2str(sum(subInds.MT))]);


%% Crop everything to 10 deg and resample MT to same elevation/abs azimuth as V1/V2

% Grab VF subset of MT within V1/V2 bounds
% get indices of MT RF's that overlap with V1/V2; allow points to be either
% left or right of fixation, at same eccentricities
% in practice, this is just the V1 range because all V2 is within V2
V1V2_x_range = max(abs([V1.x_pos V2.x_pos])); % don't care about left/right, just the max x-coordinate
V1V2_y_range = [min([V1.y_pos V2.y_pos]) max([V1.y_pos V2.y_pos])];

overlapping_RF_inds = abs(MT.x_pos) <= V1V2_x_range ...
    & MT.y_pos >= V1V2_y_range(1) ...
    & MT.y_pos <= V1V2_y_range(2);

% Cull cells outside of 10deg of eccentricities
subInds.V1 = (sqrt(V1.x_pos.^2 + V1.y_pos.^2) <= 10) & subInds.V1;
subInds.V2 = (sqrt(V2.x_pos.^2 + V2.y_pos.^2) <= 10) & subInds.V2;
subInds.MT = (sqrt(MT.x_pos.^2 + MT.y_pos.^2) <= 10) & subInds.MT & overlapping_RF_inds;

display(['V1 <= 10deg:  ' num2str(sum(subInds.V1))]);
display(['V2 <= 10deg:  ' num2str(sum(subInds.V2))]);
display(['MT <= 10deg:  ' num2str(sum(subInds.MT))]);


% Report final number of included cells
display(['Num TOTAL units after subsampling:  ' num2str(sum(subInds.V1)+sum(subInds.V2)+sum(subInds.MT))]);

keyboard