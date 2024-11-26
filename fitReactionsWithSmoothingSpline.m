function [reactionTable, data, exitCondition] = fitReactionsWithSmoothingSpline(data,optionsArray)
%fitReactions fits the reactions in the data array to the Arrenhius form
% Inputs:
% data: a struct of data that conatains the reactions. The output
% from the compareChemistryFields function or the importCrossSectionData
% function
% optionsArray: an array of booleans and initial conditions prescribed
% in the following order:
% optionsArray(1) = saveData flag
% optionsArray(2) = reactionToPlot integer
% optionsArray(3,4,5) = initalGuess array, [A, Ea, n]
% optionsArray(6) = overwriteIC flag
% optionsArray(7) = saveFMinSearch flag
% optionsArray(8) = plotIntermediate flag
% optionsArray(9) = silenceOutput flag
% optionsArray(10) = plotMultipleEEDF

% Error checking
% check is the flag options are booleans

% check if options 3-5 are numeric

% check that data is a 1x1 struct

% check that the fields of data.reactionField is either 8x1 or 3x1
% if 3x1, throw an internal flag that initial conditions are not defined
% else, use the provided inital conditions

pause(1);
%% Options Assignment
saveData = optionsArray{1};
reactionToPlot = optionsArray{2};
initialGuess = optionsArray{3};
overwriteIC = optionsArray{4};
saveFMinSearch = optionsArray{5};

plotIntermediate = optionsArray{6};
silenceOutput = optionsArray{7};
plotMultipleEEDF = optionsArray{8};
fitType = optionsArray{9};
plotTrapz = optionsArray{10};
shiftFactor = optionsArray{11};
upperLimit = optionsArray{12};

%% Declare Constants
reactionField = strcat(['Reaction_',num2str(reactionToPlot)]);
reactionName = data.(reactionField).chemistry;
% electron energy level [eV]
epsilon = data.(reactionField).data(:,1);
% cross sectional data [m2]
sigma = data.(reactionField).data(:,2);
% boltzman
kB = 8.617333262e-5; % ev/K
% mass electron
mE = 9.1093837e-31; % kilograms
% convert eV to J
eVtoJ = 1.60218e-19;


%% Declare Functions
maxTe = 70000;
minTe = 300;
maxEEDF = 2*max(epsilon);
maxPlottedEEDF = 10;
maxPlottedTe = 3000;

eedf = @(epsilon, Te) sqrt(epsilon./pi) .* 2*(1/kB/Te)^(3/2) .* exp(-epsilon/kB/Te);
epsilon_EEDF = logspace(0, log10(maxEEDF), 1000)'-1;
% Defining reaction Rate Equation
k = @(sigma, epsilon, Te) sigma .* eedf(epsilon, Te) .* sqrt(2*epsilon*eVtoJ/mE);


%% output reaction
fprintf('\n')
fprintf('Evaluating Reaction %.d: %s\n', reactionToPlot, data.(reactionField).chemistry)



%% cross sectional data
if plotIntermediate
    tile = tiledlayout;
    tileFig = tile.Parent;
    % half full screen left
    tileFig.Position = [1          53        1280         924];
    nexttile
    plotXY(epsilon, sigma, [reactionName, ': Cross Section Data'],'Electron Energy [eV]','Cross Section [m^2]')
    hold on

    [spline_Epsilon, spline_Sigma] = generateSmoothingSpline(data.(reactionField));

    plot(spline_Epsilon,spline_Sigma,'LineWidth',2)
    legend('Data','Smoothing Spline', 'Location','southeast')
    % axis padded
    % xlim([min(epsilon) max(epsilon)]);
    % ylim([min(sigma) max(sigma)]);

end



%% electron temperature evolution
load MASTER_electronTempTimeEvolution.mat time T_Te
% time is in seconds (first column), T_Te is in Kelvin (second column)
Te_evolution = [time; T_Te]';
% plot te_evolution
if plotIntermediate
    nexttile
    plotXY(time, T_Te, 'Electron Temperature Evolution Post Filament', 'Time [s]','Electron Temperature [K]')
end

% what if I choose a different Te_evolution
T_Te = linspace(minTe,maxTe,1000);
time = zeros(1,length(T_Te));
Te_evolution = [time; T_Te]';
%% EEDF
% the normalized tolerance will check if the integral of the EEDF is this
% far away from a value if 1, indiciating (if true) that is's normalized
normalizedTolerance = 1e-10;
% compute integral of EEDF to check if it's normalized
integralEEDF = integral(@(epsilon) eedf(epsilon,max(T_Te)), 0, Inf);

if (1-integralEEDF) > normalizedTolerance
    err = (1-integralEEDF);
    fprintf('Error is %.d\n',err);
    error('EEDF is not normalized');
end

if plotIntermediate
    nexttile
    plotXY(epsilon_EEDF, eedf(epsilon_EEDF,10000), [reactionName, ': EEDF for T_e = 10,000K'], 'Electron Energy [eV]','Electron Distribution [eV^-1]')
    legend([num2str(10000, '%.1f'),'K'])
    xlim([0 10])
end

if plotMultipleEEDF
    nexttile
    hold on
    plotXY(epsilon_EEDF, eedf(epsilon_EEDF,max(T_Te)), [reactionName, ': EEDF for Various T_e'], 'Electron Energy [eV]','Electron Distribution [eV^-1]')
    legend([num2str(max(T_Te), '%.1f'),'K'])
    for i = [2 5 100]
        temperature = [num2str(floor(Te_evolution(i,2))), 'K'];
        plot(epsilon_EEDF,eedf(epsilon_EEDF,Te_evolution(i,2)), '*-', 'DisplayName',temperature);
    end
    hold off
    xlim([0 10])
end

%% Reaction Rate, k
% plotting what's happening inside of k by graphing both EEDF and cross section on same axis
if plotIntermediate
    % plot product of EEDF and cross section
    nexttile
    plotXYY(epsilon_EEDF, eedf(epsilon_EEDF,10000), epsilon, sigma, 'Product of EEDF and Cross Section', 'Electron Energy [eV]','Electron Distribution [eV^-1]','Cross Section [m^2]')
    xlim([0 maxPlottedEEDF])

    hold on
    plot(spline_Epsilon,spline_Sigma,'LineWidth',2, 'Color',"#7E2F8E", 'LineStyle','-')
    ylim([0 Inf])
    legend('Normalized EEDF','Normalized Cross Section', 'Spline', 'Location','northwestoutside')


    % hand over eps/sig to the spline version
    epsilon = spline_Epsilon;
    sigma = spline_Sigma;
    % plot product of EEDF and cross section
    nexttile
    plotXY(epsilon,(sigma.*eedf(epsilon, max(T_Te))),'Product of EEDF and Cross Section','Electron Energy [eV]','[m^2 / eV]')
    ax = gca;
    ax.Children.Color = '#EDB120';
    subtitle(data.(reactionField).chemistry)
    legend('Product')
    % xlim([0 maxPlottedEEDF])
    % ylim([0 Inf])


    % plotting the k function
    if plotTrapz
        nexttile
        plot(epsilon, k(sigma, epsilon, max(T_Te)), 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b')
        title('Reaction Rate as a Function of Energy Level with Trapezoidal Integration');
        % subtitle(data.(reactionField).chemistry)
        ylabel('[m^3 * eV /s]')
        xlabel('Electron Energy [eV]')
        grid on
        xlim([0 max(epsilon)])
        hold on
        fontsize(15, "points")

        x = epsilon;
        y = k(sigma, epsilon, max(T_Te));
        % Initialize cumulative sum
        cumulative_sum = 0;

        % Plot the trapezoids used in the integration
        for i = 1:length(x)-1
            % Coordinates of the trapezoid corners
            x_trap = [x(i), x(i+1), x(i+1), x(i)];  % x points of the trapezoid
            y_trap = [0, 0, y(i+1), y(i)];          % y points of the trapezoid

            % Compute the area of the trapezoid
            trapezoid_area = (x(i+1) - x(i)) * (y(i+1) + y(i)) / 2;

            % Update cumulative sum
            cumulative_sum = cumulative_sum + trapezoid_area;

            % Fill the trapezoid with transparency
            fill(x_trap, y_trap, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 1);

            % % Calculate the position for the text (middle of the trapezoid)
            % x_text = (x(i) + x(i+1)) / 2;
            % y_text = max(y(i), y(i+1)) / 2;  % Place text midway in y-direction of the trapezoid
            % % Display the cumulative sum inside the trapezoid every 10th traps
            % if mod(i,10)==0
            %     text(x_text, y_text, sprintf('%.2e', cumulative_sum), 'FontSize', 10, 'HorizontalAlignment', 'center');
            % end
        end


        % Plot Options
        subtitle(['Trapz = ', num2str(cumulative_sum)])
        hold off;
        legend('Data', 'Trapezoids');
        fontsize(15, "points")
    end
end

% compute reaction rate
Rxn_Rate = zeros(1,size(Te_evolution,1))';

for i = 1:size(Te_evolution,1)
    Te = Te_evolution(i,2);
    Rxn_Rate(i) = trapz(epsilon, k(sigma, epsilon, Te));
    % animate reaction rate product
    % figure;
    % rxnPlot = axes;
    % ylim([0 1.5])
    % x = [epsilon,sigma];
    % sigmaFun = @(epsilon) x(ismember(x(:,1),epsilon),2);
    % k_Eps = @(epsilon) sigmaFun(epsilon) .* eedf(epsilon, Te) .* sqrt(2*epsilon*eVtoJ/mE);
    % p(1) = plot(rxnPlot,epsilon, rescale(sigma));
    % hold on
    % p(2) = plot(rxnPlot,epsilon, rescale(eedf(epsilon,Te)));
    % p(3) = plot(rxnPlot,epsilon, rescale(k(sigma, epsilon, Te)));
    % p(4) = plot(rxnPlot,epsilon, rescale(k_Eps(epsilon)));
    % hold off
    % xlabel('Epsilon, [eV]')
    % title(['Te: ',num2str(Te, '%.2f'),'K'])
    % legend('Cross Section','EEDF','K')
    % pause(0.1);
    % delete(p)
    % drawnow

end


% convert RxnRate from m3/s to cm3/s
Rxn_Rate = Rxn_Rate .* 100^3;

if min(Rxn_Rate) < 0
    disp('Negative Reaction Rate Detected')
    reactionTable = manualReactionTable(data);
    exitCondition = 'negReactionRate';
    return
end

if plotIntermediate
    nexttile;
    plotLogLog(T_Te, Rxn_Rate, [reactionName, ': Reaction Rates'], 'Electron Temperature [K]','Reaction Rate [cm^3 / s]' )
    ax = gca;
    ax.Children.LineStyle = 'none';
end


%% Curve Fitting
% Define the objective function to minimize
% transpose to make like rest of variables
% fitting_Te = Te_evolution(:,2);
% I want to prioritize fitting in the lower regieme of temperatures
% lowerLimit = 80;
% if isfield(data.(reactionField), 'lowerLimit')
%     upperLimit = data.(reactionField).lowerLimit;
% else
% end


fitting_Te = Te_evolution(1:upperLimit,2);
fitting_Rxn_Rate = Rxn_Rate(1:upperLimit);

objectiveFunction = @(params) computeSSE(params, fitting_Te, fitting_Rxn_Rate);

% Initial guesses for A, Ea, and n
% check if initial conditions are saved in the stucture, and if so, use
% those instead
if isfield(data.(reactionField), 'InitialConditions') && ~overwriteIC
    initialParams = data.(reactionField).InitialConditions;

    if ~silenceOutput
        disp('Loading Initial Conditions from Structure')
        disp(['Initial A: ', num2str(initialParams(1))]);
        disp(['Initial Ea: ',num2str(initialParams(2))]);
        disp(['Initial n: ', num2str(initialParams(3))]);
    end

else
    initialParams = initialGuess;  % [A, Ea, n] adjust these as needed
    if ~silenceOutput
        disp('Loading Initial Conditions from Initial Guess')
        disp(['Initial A: ', num2str(initialParams(1))]);
        disp(['Initial Ea: ',num2str(initialParams(2))]);
        disp(['Initial n: ', num2str(initialParams(3))]);
    end
end

% Set optimization options, increasing MaxFunEvals and MaxIter
n = 1e4;
options = optimset('MaxFunEvals', n, 'MaxIter', n, 'Display','off');  % Increase as needed

% Use fminsearch to minimize the SSE
optimizedParams = fminsearch(objectiveFunction, initialParams, options);
% Display the optimized parameters
if ~silenceOutput
    A_opt = optimizedParams(1);
    Ea_opt = optimizedParams(2);
    n_opt = optimizedParams(3);
    disp(['Optimized A: ', num2str(A_opt)]);
    disp(['Optimized Ea: ', num2str(Ea_opt)]);
    disp(['Optimized n: ', num2str(n_opt)]);
end
final_sse = computeSSE(optimizedParams, fitting_Te, fitting_Rxn_Rate);


modelfun = @(b,x) b(1).*x(:,1).^b(3) .* exp(-b(2)./x(:,1));
Plot_RxnRateComp = figure;
% half full screen right
Plot_RxnRateComp.Position = [1282          53        1279         924];
% true data points
% semilogy(T_Te, Rxn_Rate, '*')
plot(T_Te, Rxn_Rate, '*')
hold on
modelTe = Te_evolution(:,2)';
% modelTe = 3000:10:60000;
% plot fminsearch
plot(modelTe,modelfun(optimizedParams, modelTe'), '--', 'LineWidth',2)

Plot_RxnRateComp.CurrentAxes.Title.String = strcat([data.(reactionField).chemistry,' Fit Comparison']);
Plot_RxnRateComp.CurrentAxes.XLabel.String = 'Electron Temperature [K]';
Plot_RxnRateComp.CurrentAxes.YLabel.String = 'Reaction Rate [cm^3 / s]';
Plot_RxnRateComp.CurrentAxes.XAxis.Exponent = 0;

xline(Te_evolution(upperLimit,2), 'LineWidth',2)
legend('Data', 'Nonlinear Fitting','Fitting Limit','Location','southeast')
fontsize(15, "points")
grid on
xlim([minTe Te_evolution(upperLimit,2)*1.5])
hold off

%% Save Data

if saveData
    data.(reactionField).InitialConditions = initialParams;

    data.(reactionField).A = A_opt;
    data.(reactionField).Ea = Ea_opt;
    data.(reactionField).n = n_opt;
    data.(reactionField).SSE = final_sse;

    data.(reactionField).shiftFactor = shiftFactor;
    data.(reactionField).crossSectionFitTypes = fitType;

    data.(reactionField).lowerLimit = upperLimit;


    % save('BOLSIG_Rates_Data',"data")
end

%% Create Reaction Table
% Assume 'data' is the input 1x1 structure with 67 fields.
% Each field has subfields '.chemistry', '.A', '.Ea', and '.n'.

% Get all field names of the structure
fields = fieldnames(data);

% Preallocate cell arrays for the table columns
fieldNames = cell(length(fields), 1);
chemistryValues = cell(length(fields), 1);
A_values = zeros(length(fields), 1);  % Assuming numeric values for A
Ea_values = zeros(length(fields), 1); % Assuming numeric values for Ea
n_values = zeros(length(fields), 1);  % Assuming numeric values for n

% Loop over all fields and extract the relevant values
for i = 1:length(fields)
    fieldName = fields{i};
    fieldNames{i} = fieldName;  % Store the field name
    % if the fitting exists, that is, if the data.(fieldName) length is
    % greater than 3, save the data, else, input NaN
    if length(fieldnames(data.(fieldName))) > 3
        chemistryValues{i} = data.(fieldName).chemistry;  % Store the chemistry value
        A_values(i) = data.(fieldName).A;  % Store the A value
        Ea_values(i) = data.(fieldName).Ea;  % Store the Ea value
        n_values(i) = data.(fieldName).n;  % Store the n value
    else
        chemistryValues{i} = data.(fieldName).chemistry;  % Store the chemistry value
        A_values(i) = NaN;  % Store the A value as NaN if no fit exists
        Ea_values(i) = NaN;  % Store the Ea value as NaN if no fit exists
        n_values(i) = NaN;  % Store the n value as NaN if no fit exists
    end
end

% Create a table with five columns: 'ReactionName', 'Chemistry', 'A', 'Ea', and 'n'
reactionTable = table(fieldNames, chemistryValues, A_values, Ea_values, n_values, ...
    'VariableNames', {'ReactionName', 'Chemistry', 'A', 'Ea', 'n'});

% Display the table
% disp(reactionTable);
% save("reactionTable",'reactionTable')

%     proceed = input('Press Enter to proceed, or any other key to exit: ', 's');
%     if isempty(proceed)  % If Enter is pressed (empty input)
%         % disp(['Continuing the loop, iteration: ', num2str(i)]);
%         % % Your code for this iteration
%         delete(Plot_RxnRateComp)
%         drawnow
%     else
%         % disp('Exiting the loop.');
%         break;
%     end
%
% end
disp(data.(reactionField));



% end of function
exitCondition = 'normal';
end

function sse = computeSSE(params, fitting_Te, Rxn_Rate)
A = params(1);
Ea = params(2);
n = params(3);

% Arrhenius equation
k_fit = A .* fitting_Te.^n .* exp(-Ea ./ fitting_Te);

% Calculate the sum of squares due to error (SSE)
sse = sum((Rxn_Rate - k_fit).^2);
end