function [epsilonFit,sigmaFit] = generateSmoothingSpline(inputStruct)
% smoothing spline fit
debugMode = false;
% get data extracted from input
XC = inputStruct.data;
epsilon = XC(:,1);
sigma = XC(:,2);
maxEps = max(epsilon);
minEps = min(epsilon);
% if minEps is already 0, get the next biggest value from
while(min(minEps)<=0)
    if minEps == 0
        minEps = mink(epsilon,2);
        minEps = max(minEps);
    end
end

% check for trailing 0's in epsilon, ID first one
trail = sigma ~=0;
index = find(trail ~= 0, 1, 'last');
if index == length(epsilon)
    maxEps = epsilon(index);
else
    maxEps = epsilon(index+1);
end

% check for leading 0's in epsilon, ID first one
trail = sigma ~=0;
index = find(trail ~= 0, 1, 'first');
if index == 1
    minEps = epsilon(index);
else
    minEps = epsilon(index-1);
end


% plotting
if debugMode
    figure;
    plot(epsilon,sigma,'*');
    hold on
end


% remove zeros from sigma
r = sigma>0;
sigma = sigma(r);
epsilon = epsilon(r);

% remove left edge outliers
[maxSig, idxMaxSig] = max(sigma);
[maxSig2, idxMaxSig2] = max(sigma(idxMaxSig+1:end));
if idxMaxSig2 ~= 1
    sigma(idxMaxSig) = min(sigma);
    warning('Found an outlier in sigma data. Check fit')
end

if debugMode
    plot(epsilon,sigma,'*');
    hold on
end


% define a fitting version of the input variables
% inside of log space
fittingX = log(epsilon);
fittingY = -log(sigma);




% function to generate smoothing spline with a smoothing parameter of 1
[fitResult,~] = createFit(fittingX, fittingY);

% define the fitting regime based on inputs
if minEps <= 0
    epsilonFit = log(minEps:0.1:maxEps*1.2)';
else
    epsilonFit = log(0.8*minEps:0.1:maxEps*1.2)';
end

% get fitted result in log space
sigmaFit = fitResult(epsilonFit);

% convert fitted result back into normal space
sigmaFit = exp(-sigmaFit);
epsilonFit = exp(epsilonFit);


% find the right edge of the fitting data
rightEdgeValue = max(exp(fittingX));
x = find(epsilonFit>=rightEdgeValue);
rightEdgeFit = x(1);
% scale to the fit data

% [midpoint] = round(length(sigmaFit)/2);

% find second IP
% loop through diffs
diffs = diff(sigmaFit(rightEdgeFit:end));
foundIP = false;
for i = 2:length(diffs)
    % if the d/dx is positive, and all the previous values was negative
    % it's the inflection point
    if diffs(i) > 0 && diffs(i-1) <= 0
        IP = i+rightEdgeFit;
        sigmaFit(IP:end) = 2*sigmaFit(IP) - sigmaFit(IP:end);
        foundIP = true;
        break;
    end
end
% if it goes negative, force it to 0
sigmaFit(sigmaFit<0) = 0;

if debugMode
    xlim([0 50])
    plot(epsilonFit(IP),sigmaFit(IP), '*');
    plot(epsilonFit,sigmaFit, 'LineWidth',2);
end

% find the left edge of the fitting data
leftEdgeValue = min(exp(fittingX));
x = find(epsilonFit<=leftEdgeValue);
leftEdgeFit = x(1);

% find first IP
diffs = diff(sigmaFit(1:leftEdgeFit));

for i = 2:length(diffs)
    if diffs(i) > 0 && diffs(i-1) <= 0
        IP = i;
        % if we already found an IP, apply transform to first half of data
        % else, apply to second half
        if foundIP
            % flip Y data over a YLILE defined at the IP
            sigmaFit(1:IP) = 2*sigmaFit(IP) - sigmaFit(1:IP);
        else
            sigmaFit(IP:end) = 2*sigmaFit(IP) - sigmaFit(IP:end);
        end
        break;
    end
end

if debugMode
    plot(epsilonFit(IP),sigmaFit(IP), '*');
end
% if it goes negative, force it to 0
sigmaFit(sigmaFit<0) = 0;


if debugMode
% xlim([0 50])
    plot(epsilonFit,sigmaFit, 'LineWidth',2);
end
% force to 0 after max and min Eps
sigmaFit(epsilonFit<=minEps) = 0;
sigmaFit(epsilonFit>=maxEps) = 0;

% sigmaFit(epsilonFit<=min(epsilon)) = 0;

if debugMode
    % plotting
    plot(epsilonFit,sigmaFit, 'LineWidth',2);
    title(inputStruct.chemistry)
    xlabel( 'epsilon', 'Interpreter', 'none' );
    ylabel( 'sigma', 'Interpreter', 'none' );
    grid on
    legend('Experimental Data','Curve Fit')
    % % ylim([0 max(sigma)])
end



end