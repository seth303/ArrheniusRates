%% Import Table
clc; clear; %close all; 
% define the git repo
repo = gitrepo('ZDPlasKin FS Simulations/');
currBranch = repo.CurrentBranch.Name;
% define filepaths
arrRatesFilepath =     'ZDPlasKin FS Simulations/Results/NRP_300um_1200K_FIELD0.02/qt_rates.txt';
arrReactionListFilepath = 'ZDPlasKin FS Simulations/Results/NRP_300um_1200K_FIELD0.02/qt_reactions_list.txt';
bolsigReactionListFilepath = 'ZDPlasKin FS Simulations/qt_reactions_list_bolsig.txt';
temperatureFilePath =  'ZDPlasKin FS Simulations/Te_Arrenhius';
bolsigRatesFilePath =  'ZDPlasKin FS Simulations/qt_rates_bolsig+.txt';
% ID import options, and set them
opts = detectImportOptions(arrRatesFilepath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
% import arrhenius rates
arrheniusRates = readtable(arrRatesFilepath,opts);
% import temperature
opts = detectImportOptions(temperatureFilePath);
electronTemp = readtable(temperatureFilePath,opts);
electronTemp.Properties.VariableNames = {'Time [s]','Electron Temp [K]'};
% import reaction names
varNames = readcell(arrReactionListFilepath);
varNames = [{0,'Time [s]'};varNames];
varNames = addSuffixToDuplicates(varNames);
% combine rates and reaction names
arrheniusRates.Properties.VariableNames = varNames(:,2);
% % assign the table to the git repo structure
% data.(currBranch) = T;
% extract the table (with the git repo as the name)
% struct2vars(data);
% clear the old table and placeholder structure
clear T data
% import bolsig rates
opts = detectImportOptions(bolsigRatesFilePath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
varNames = readcell(bolsigReactionListFilepath);
varNames = [{0,'Time [s]'};varNames];
BOLSIGRates = readtable(bolsigRatesFilePath,opts);
BOLSIGRates.Properties.VariableNames = varNames(:,2);
% data.BOLSIGRates = T;
% struct2vars(data);
clear T data
% load master reaction table
load ../MASTER_reactionTable.mat
load ../MASTER_BOLSIG_Rates_Data.mat
% import species densitites
% define filepaths
arrDensitiesFilePath =     'ZDPlasKin FS Simulations/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt';
speciesListFilepath =      'ZDPlasKin FS Simulations/Results/NRP_300um_1200K_FIELD0.02/qt_species_list.txt';
bolsigDensitiesFilePath =  '../Isolate Reactions/main/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt';


% ID import options, and set them
opts = detectImportOptions(arrDensitiesFilePath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
% import arrhenius densities
arrheniusDensities = readtable(arrDensitiesFilePath,opts);

% import reaction names
varNames = readcell(speciesListFilepath);
varNames = [{0,'Time [s]'};varNames];
% combine rates and reaction names
arrheniusDensities.Properties.VariableNames = varNames(:,2);
% assign the table to the git repo structure
% data.(currBranch) = T;
% % extract the table (with the git repo as the name)
% struct2vars(data);
% clear the old table and placeholder structure

% import bolsig rates
opts = detectImportOptions(bolsigDensitiesFilePath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
BOLSIGDensities = readtable(bolsigDensitiesFilePath,opts);
BOLSIGDensities.Properties.VariableNames = varNames(:,2);
% data.BOLSIG_Densities = T;
% struct2vars(data);
%% Reaction Number
clc; close all;
plotCombined = 0;
scale = 1e9;
misMatches = zeros(size(arrheniusRates,2),1);
for rxn = 2:size(arrheniusRates,2)

    % arrhenius output
    k_arr = table2array(arrheniusRates(:,rxn));

    % BOLSIG+
    n_bolsig = table2array(BOLSIGRates(:,rxn));

    % find differences in k
    tol = 0.999289696677690*0.9;
    correlation = corr(n_bolsig, k_arr);
    if correlation > tol
        misMatches(rxn) = rxn-1;
    end
end
% find bolsig reactions
bolsigIdx = contains(string((varNames(:,2))), 'bolsig');
rxnNames = string((varNames));
bolsigRxnsNum = rxnNames(bolsigIdx);
bolsigRxnsName = rxnNames(str2double(bolsigRxnsNum)+1,2);

% if plotCombined
%     tL = tiledlayout;
% else
%     fig = figure;
% end
Te_arr = table2array(electronTemp(:,2));
[M,I] = max(Te_arr);
% N2->N2(C3)
% bolsigRxnsNum = "74";
% bolsigDensityNum = "14";

bolsigRxnsNum = ["65","66","67", "68","69","70","74","75","118"]';
bolsigDensityNum = ["11","11","11","12","12","12","14","14","39"];
electronDensityNum = 56;
iter = 1;
for rxn = (str2double(bolsigRxnsNum))'
    % rxn = 65;
    if plotCombined
        nexttile
    else
        fig = figure;
    end
    % plot arrhenius output (1/2) reaction rate
    % arrheniusRates is a table with rates with units of cm^-3 s^-1
    n_arr =  table2array(arrheniusRates(:,rxn+1));
    disp(['Found Arrhenius Rate [1/(cm^3*s^3)] for reaction: ',char(arrheniusRates(:,rxn+1).Properties.VariableNames)])
    % plot(Te_arr(I:end), n_arr(I:end), '*');
    hold on
    % convert reaction rate to reaction rate constant
    % reactionRate[cm^-3 s^-1]
    % reactionRateConstant [cm^3/s] so we need to multiply by 
    % dn_a/dt = k *n_a*n_elec
    % k = dn_a/dt/n_a/n_elec
    speciesDens_arr = table2array(arrheniusDensities(:,str2double(bolsigDensityNum(iter))));
    disp(['Found Arr. Species Density [1/(cm^3)] for species: ',char(arrheniusDensities(:,str2double(bolsigDensityNum(iter))).Properties.VariableNames)])
    elecDens_arr = table2array(arrheniusDensities(:,electronDensityNum));
    disp(['Found Arr. Species Density [1/(cm^3)] for species: ',char(arrheniusDensities(:,electronDensityNum).Properties.VariableNames)])
    % cm^3/s = (1/cm^3 *s)/(1/cm^3)/(1/cm^3)
    k_arr = n_arr./speciesDens_arr./elecDens_arr;
    plot(Te_arr(I:end), k_arr(I:end), 'o');
    % plot(Te_arr, k_arr, 'o');


    % plot BOLSIG+ (1/2)
    n_bolsig = table2array(BOLSIGRates(:,rxn+1));
    disp(['Found BOLSIG Rate [1/(cm^3*s^3)] for reaction: ',char(BOLSIGRates(:,rxn+1).Properties.VariableNames)])

    % plot(Te_arr(I:end), n_bolsig(I:end), '^');

    % convert reaction rate to reaction rate constant
    speciesDens_bolsig = table2array(BOLSIGDensities(:,str2double(bolsigDensityNum(iter))));
    disp(['Found BOLISG Species Density [1/(cm^3)] for species: ',char(BOLSIGDensities(:,str2double(bolsigDensityNum(iter))).Properties.VariableNames)])
    elecDens_bolsig = table2array(BOLSIGDensities(:,electronDensityNum));
    disp(['Found BOLSIG Species Density [1/(cm^3)] for species: ',char(BOLSIGDensities(:,electronDensityNum).Properties.VariableNames)])

    k_ct_bolsig = n_bolsig./speciesDens_bolsig./elecDens_bolsig;
    plot(Te_arr(I:end), k_ct_bolsig(I:end), '+');


    % seach for matching reaction in reactionTable
    reactionName_Arrhenius = upper(string(BOLSIGRates.Properties.VariableNames(rxn+1)));
    % reformat this string to match the format of reactionTable
    reactionName_Arrhenius = strrep(reactionName_Arrhenius,'=>','->');
    if contains(reactionName_Arrhenius, 'BOLSIG:')
        title(reactionName_Arrhenius)
        reactionName_Arrhenius = strsplit(reactionName_Arrhenius, 'BOLSIG:');
        reactionName_Arrhenius = reactionName_Arrhenius(2);
    else
        title(reactionName_Arrhenius)
    end

    disp(strcat("Looking for: ", reactionName_Arrhenius))
    % searchString = strsplit(reactionName_BOLSIG, 'bolsig:');
    % searchString = searchString{2};
    % % remove single '
    % searchString = strrep(searchString, '''', '`');
    % remove white space
    searchArray  = upper(regexprep(reactionTable.Chemistry, '\s+', ''));
    % lowercase all letters in the parenthesis
    % searchString = regexprep(searchString, '\(([^)]+)\)', '${lower($0)}');
    % searchString = 'O2(a1D)->O2^+'
    [TF, Idx] = ismember(searchArray, reactionName_Arrhenius);
    arrheniusReaction = reactionTable(TF,:);

    if isempty(arrheniusReaction)
        % % try the backwards arrow
        % disp('None found. Looking for alternative versions')
        % disp('Checking <->')
        % reactionName_Arrhenius = strrep(reactionName_Arrhenius,'->','<->');
        % [TF, Idx] = ismember(searchArray, reactionName_Arrhenius);
        % arrheniusReaction = reactionTable(TF,:)
        warning('Cannot find reaction in Reaction Table. Check search string')
    else
        % plot arrhenius form
        A = arrheniusReaction.A;
        Ea = arrheniusReaction.Ea;
        n = arrheniusReaction.n;
        arrheniusForm = @(Te) scale*A.*Te.^n.*exp(-Ea./Te);
        plot(Te_arr(I):-1:Te_arr(end), arrheniusForm(Te_arr(I):-1:Te_arr(end)))
    end


    % % fminsearch for BOLSIG data
    % % Set optimization options, increasing MaxFunEvals and MaxIter
    % n = 1e4;
    % options = optimset('MaxFunEvals', n, 'MaxIter', n, 'Display','off','PlotFcns',[]);  % Increase as needed
    % lowerLimit = 528;
    % fitting_Te = Te_arr(lowerLimit:end);
    % fitting_Rxn_Rate = k_ct_bolsig(lowerLimit:end);
    % fitting_Rxn_Rate(fitting_Rxn_Rate==Inf)=0;
    % initialParams = data.(string(arrheniusReaction.ReactionName)).InitialConditions;
    % objectiveFunction = @(params) computeSSE(params, fitting_Te, fitting_Rxn_Rate);
    % 
    % % Use fminsearch to minimize the SSE
    % 
    % optimizedParams = fminsearch(objectiveFunction, initialParams, options);
    % 
    % modelfun = @(b,x) b(1).*x.^b(3) .* exp(-b(2)./x);
    % modelTe = 0:1:Te_arr(I);
    % plot(modelTe,modelfun(optimizedParams, modelTe'), 'LineWidth',2)

    if ~plotCombined
        % plot options
        xlabel('Electron Temp [K]')
        ylabel('Reaction Rate [cm^3/s]')
        fig.Children.YScale = 'log';
        fig.Children.XAxis.Exponent = 0;
        fontsize(20, "points")
        legend('Arrhenius Rate Ct.','BOLSIG+ Rate Ct.', 'Arrhenius Curve Fit', 'FMinSearch', 'Location','best')
        % xlim([0 61000])
    end
    iter = iter+1;
    hold off
    grid on
end



% plot options
if plotCombined
    legend('Arrhenius Rate','Arrhenius Rate Ct.','BOLSIG+ Rate','BOLSIG+ Rate Ct.', 'Arrhenius Curve Fit', 'FMinSearch')
    xlabel(tL, 'Electron Temp [K]', 'FontSize',15)
    ylabel(tL, 'Reaction Rate [cm^3/s]', 'FontSize',15)
    for i = 2:size(tL.Children,1)
        tL.Children(i).YScale = 'log';
        tL.Children(i).XAxis.Exponent = 0;
    end
    tL.Children(1).FontSize = 14;
    tL.Children(1).Position = [0.837252869897959,0.106846473029046,0.108617665816327,0.119232365145228];
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