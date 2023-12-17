function [retDisps] = screen2retDisp(scrDisps,rfCent,cortArea)

% Convert frontoparallel "screen" disparity to actual retinal disparity

switch cortArea
    case 'V1'
        z   = 0.89;
    case 'V2'
        z   = 0.89;
    case 'MT'
        z   = 0.57;
end

% assume a 6cm ipd
ipd = 0.03;
xL  = -ipd/2;
xR  = ipd/2;

retDisps = zeros(1,numel(scrDisps));

% Define equations
posL   = @(z,xL,xP) atan2d(xP-xL,z) - atan2d(xL,z);
posR   = @(z,xR,xP) atan2d(xP-xR,z) - atan2d(xR,z);
dotSep = @(z,disps,xR) 2*(z*tand(disps/2 + atan2d(xR,z)) - xR);
scrPos = @(z,angs) z*tand(angs);

% Calculate true disparity
for ii = 1:numel(scrDisps)

    % Get dot positions on screen
    xPos = scrPos(z,rfCent) + dotSep(z,scrDisps(ii),xR)/2*[-1;1];

    % Get eye positions
    eyeL = posL(z,xL,xPos(1));
    eyeR = posR(z,xR,xPos(2));

    % Get retinal disparities
    retDisps(ii) = eyeR - eyeL + 2*atan2d(xR,z) - 2*atan2d(xL,z);

end
