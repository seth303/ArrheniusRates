% smoothing spline fit
clc; close all; clear
load MASTER_BOLSIG_Rates_Data.mat

% get data
XC = data.Reaction_157.data;
data.Reaction_157.chemistry
epsilon = XC(:,1);
sigma = XC(:,2);
% remove zeros from sigma
r = sigma>0;
sigma = sigma(r);
epsilon = epsilon(r);

fittingX = log(epsilon);
fittingY = -log(sigma);
fittingY(fittingY==Inf) = 0;

[fitResult,gof] = createFit(fittingX, fittingY);

% Plot fit with data.
figure( 'Name', 'N2 -> N2(A3,v0-4) Original Data' );
h = plot(epsilon,sigma,'*');
% legend( h, 'fittingY vs. fittingX', 'N2 -> N2(A3,v0-4)', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'epsilon', 'Interpreter', 'none' );
ylabel( 'sigma', 'Interpreter', 'none' );
grid on
hold on
if max(XC(:,1)) < 500
epsilon = log(0:0.1:max(XC(:,1)))';
else
    epsilon = log(0:0.1:200)';
end
sigmaFit = fitResult(epsilon);


% sigmaFit(sigmaFit<min(fittingY)) = 0;
% sigmaFit(sigmaFit>max(fittingY)) = 0;

sigmaFit = exp(-sigmaFit);

% if the derivative of sigma after the last fitting point turns positive,
% force it to 0
% find the inflection point in the second half of the data
[I] = round(length(sigmaFit)/2);
diffs = diff(sigmaFit(I:end));

% loop through diffs
for i = 2:length(diffs)
    % if the d/dx is positive, and all the previous values was negative
    % it's the inflection point
    if diffs(i) > 0 && diffs(i-1) <= 0
        IP = i+I;
        sigmaFit(IP:end) = 2*sigmaFit(IP) - sigmaFit(IP:end);
        plot(exp(epsilon(IP)), sigmaFit(IP),'*')
        break;
    end
end




% sigmaFit(sigmaFit == 1) = 0;
% plot(exp(epsilon),sigmaFit, 'LineWidth',2);
hold on

% from the infletion point onwards, apply a negative diffs
% sigmaFit(IP+1:end) = sigmaFit(IP+1:end) .* -diffs(IP:end);
% flip Y data over a YLILE defined at the IP

% plot(exp(epsilon),sigmaFit, 'LineWidth',2);
% if it goes negative, force it to 0
sigmaFit(sigmaFit<0) = 0;
% plot(exp(epsilon),sigmaFit, 'LineWidth',2);

% find first IP
% find the inflection point in the second half of the data
% loop through diffs
diffs = diff(sigmaFit(1:I));

for i = 2:length(diffs)
    % if the d/dx is positive, and all the previous values was negative
    % it's the inflection point
    if diffs(i) > 0 && diffs(i-1) <= 0
        IP = i;
        % flip Y data over a YLILE defined at the IP
        sigmaFit(1:IP) = 2*sigmaFit(IP) - sigmaFit(1:IP);
        break;
    end
end
plot(exp(epsilon(IP)), sigmaFit(IP),'*')



% plot(exp(epsilon),sigmaFit, 'LineWidth',2);
% if it goes negative, force it to 0
sigmaFit(sigmaFit<0) = 0;

epsilon = exp(epsilon);
plot((epsilon),sigmaFit, 'LineWidth',2);
