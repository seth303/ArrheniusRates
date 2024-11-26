%% Import Table
clc; clear; %close all;
% define the git repo
mainBranch = gitrepo('main/');
% arrheniusRatesBranch = gitrepo('arrheniusRates/');
% arrheniusRatesBranch = gitrepo('fsPostFilament/');
arrheniusRatesBranch = gitrepo('Ivanov/');


% pull(mainBranch)
% pull(arrheniusRatesBranch)
arrheniusRatesBranch.LastCommit.Message


% currBranch = arrheniusRatesBranch.CurrentBranch.Name;
% define filepaths
if contains(arrheniusRatesBranch.WorkingFolder,'arrheniusRates')
    arrDensitiesFilePath =     strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt');
    speciesListFilepath =      strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_species_list.txt');
    bolsigDensitiesFilePath =  strcat(mainBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt');
    timescale = 'ns';
end
% fs
if contains(arrheniusRatesBranch.WorkingFolder,'fsPostFilament')
    arrDensitiesFilePath =     strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_densities.txt');
    speciesListFilepath =      strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_species_list.txt');
    bolsigDensitiesFilePath =  strcat(mainBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt');
    arrRatesFilepath =         strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_rates.txt');
    arrReactionListFilepath =  strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_reactions_list.txt');
    arrConditionsFilePath =    strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_conditions.txt');
    arrConditionsListFilePath =    strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_conditions_list.txt');


    timescale = 'fs';
end

if contains(arrheniusRatesBranch.WorkingFolder,'Ivanov')
    
    arrDensitiesFilePath =     strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_densities.txt');
    speciesListFilepath =      strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_species_list.txt');
    bolsigDensitiesFilePath =  strcat(mainBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt');
    arrRatesFilepath =         strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_rates.txt');
    arrReactionListFilepath =  strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_reactions_list.txt');
    arrConditionsFilePath =    strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_conditions.txt');
    arrConditionsListFilePath =    strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_conditions_list.txt');

    timescale = 'fs';
end

% ID import options, and set them
opts = detectImportOptions(arrConditionsFilePath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
% import arrhenius rates
arrheniusConditions = readtable(arrConditionsFilePath,opts);

% import condititons names
arrheniusConditions.Properties.VariableNames(2) = {'Reduced field [Td]'};
arrheniusConditions.Properties.VariableNames(3) = {'Gas temperature [K]'};
arrheniusConditions.Properties.VariableNames(4) = {'Electron temperature [K]'};

% ID import options, and set them
opts = detectImportOptions(arrRatesFilepath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
% import arrhenius rates
arrheniusRates = readtable(arrRatesFilepath,opts);

% import reaction names
varNames = readcell(arrReactionListFilepath);
varNames = [{0,'Time [s]'};varNames];
varNames = addSuffixToDuplicates(varNames);
% combine rates and reaction names
arrheniusRates.Properties.VariableNames = varNames(:,2);

% ID import options, and set them
opts = detectImportOptions(arrDensitiesFilePath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
% import arrhenius densities
arrhenius_Densities= readtable(arrDensitiesFilePath,opts);

% import reaction names
varNames = readcell(speciesListFilepath);
varNames = [{0,'Time [s]'};varNames];
% combine rates and reaction names
arrhenius_Densities.Properties.VariableNames = varNames(:,2);


% import bolsig rates
opts = detectImportOptions(bolsigDensitiesFilePath);
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
BOLSIG_Densitites = readtable(bolsigDensitiesFilePath,opts);
BOLSIG_Densitites.Properties.VariableNames = varNames(:,2);
% load master reaction table
load ../MASTER_reactionTable.mat

clear opts varNames
%% Compare Differences
if strcmp(timescale,'ns')
    clc;
    errorCounter = 0;
    startLooking =  1;

    for n = 1:size(arrhenius_Densities,2)

        % find rows that are different:
        column_Different = find(~all(BOLSIG_Densitites{startLooking:end,n} == arrhenius_Densities{startLooking:end,n},2));

        % check if difference is less than some threshold
        err = zeros(1,length(column_Different));
        for i = 1:length(column_Different)
            expected = table2array(BOLSIG_Densitites(column_Different(i),n));
            actual = table2array(arrhenius_Densities(column_Different(i),n));

            err(i) = abs((actual-expected)/expected)*100;
        end
        [M,I] = max(err);

        if M > 10
            errorCounter  = errorCounter + 1;
            errorReactions(errorCounter) = n;
            errorLocations(errorCounter) = column_Different(I)+startLooking;
        end

    end
end
%% Plot Differences
if strcmp(timescale,'ns')
    timeBOLSIG = table2array(BOLSIG_Densitites(:,1));
    timeArrhenius = table2array(arrhenius_Densities(:,1));
    close all;
    figure;
    TL = tiledlayout;
    errorReactions = 1:1:56;
    % Nitrogen Reactions: 1-11
    % Oxygen Reactions: 12-27
    % Combo Reactions: 28-41
    maxPlottedReactions = 11;
    for i = 1:length(errorReactions)
        nexttile;

        BOLSIG = table2array(BOLSIG_Densitites(:,errorReactions(i)));
        Individual = table2array(arrhenius_Densities(:,errorReactions(i)));
        %
        % plot(localTime,BOLSIG(localIndex),'*')
        % hold on
        % plot(localTime,Individual(localIndex), '^')

        plot(timeBOLSIG,BOLSIG,'LineWidth',1)
        hold on
        plot(timeArrhenius,Individual,'LineWidth',2, 'LineStyle','--')
        title(BOLSIG_Densitites.Properties.VariableNames(errorReactions(i)))
        grid on
        if i == maxPlottedReactions
            break
        end

    end

    % lgnd = legend('BOLSIG+ Error Region','Wilson Error Region', 'BOLSIG+ Global','Wilson Global');
    lgnd = legend('BOLSIG+','Wilson');
    lgnd.Layout.Tile = 'north';
    xlabel(TL,'Time [s]')
    ylabel(TL,'Density [cm^-3]')
    fontsize(20, "points")




    for i = 1:length(TL.Children)
        if isgraphics(TL.Children(i),'Axes')
            % TL.Children(i).YScale = 'log';
        end
    end
    % %
    % TL.Parent.WindowState = 'maximized';
    % pause(0.5)
    % saveas(TL.Parent,arrheniusRatesBranch.LastCommit.Message,'png')
    % TL.Parent.WindowState = 'normal';

end
%% Plotting Individual Fs Reactions
if strcmp(timescale,'fs')
    timeArrhenius = table2array(arrhenius_Densities(:,1));
    close all;
    fsPlot = figure;
    errorReactions = 1:1:56;
    % N2(A)
    maxPlottedReactions = 999;
    for i = [11 14 21 2 56]
        Individual = table2array(arrhenius_Densities(:,errorReactions(i)));

        hold on
        plot(timeArrhenius,Individual,'LineWidth',2, 'LineStyle','-')
        % title(BOLSIG_Densitites.Properties.VariableNames(errorReactions(i)))
        grid on
        if i == maxPlottedReactions
            break
        end

    end

    % lgnd = legend('BOLSIG+ Error Region','Wilson Error Region', 'BOLSIG+ Global','Wilson Global');
    lgnd = legend('N2(A)','N2(C3)','N4^+','N2','e-');
    xlabel('Time [s]')
    ylabel('Density [cm^-3]')
    fontsize(20, "points")
    fsPlot.Children(2).YScale = 'log';
    xlim([0 10e-9])

    % %
    % TL.Parent.WindowState = 'maximized';
    % pause(0.5)
    % saveas(TL.Parent,arrheniusRatesBranch.LastCommit.Message,'png')
    % TL.Parent.WindowState = 'normal';

end

%% Plotting indivdiual reaction rates
% close all
if strcmp(timescale,'fs')
    RxnsNum = ["3","77","197"]';
    % DensityNum = ["21","11","11","12","12","12","14","14","39"];
    electronDensityNum = 56;
    iter = 1;
    Te_arr = table2array(arrheniusConditions(:,4));
    [M,I] = max(Te_arr);
    time = table2array(arrheniusRates(:,1));
    fig = figure;
    hold on
    for rxn = (str2double(RxnsNum))'
        % plot arrhenius output (1/2) reaction rate
        % arrheniusRates is a table with rates with units of cm^-3 s^-1
        n_arr =  table2array(arrheniusRates(:,rxn+1));
        disp(['Found Arrhenius Rate [1/(cm^3*s^3)] for reaction: ',char(arrheniusRates(:,rxn+1).Properties.VariableNames)])

        % if rxn is the N2(E3), add them together
        if rxn == 77
            n_arr = n_arr +  table2array(arrheniusRates(:,78+1));
        end
        % k_arr = n_arr./speciesDens_arr./elecDens_arr;
        plotTitle = char(arrheniusRates(:,rxn+1).Properties.VariableNames);
        % remove the 3 from C3
        plotTitle = erase(plotTitle, '3');
        plot(time, n_arr, 'LineWidth',2, 'DisplayName',plotTitle);
    end
    hold off
    xlim([0 4e-9])
    fig.Children.YAxis.Scale = 'log';
    fontsize(15,'points')
    xlabel('Time [s]')
    ylabel('Reaction Rate [1/cm^3*s]')
    legend
    title('N2(C) Production Mechanism')
    grid on
     % N2(C) losses
    RxnsNum = ["202","203","173"]';
    % DensityNum = ["21","11","11","12","12","12","14","14","39"];
    electronDensityNum = 56;
    iter = 1;
    Te_arr = table2array(arrheniusConditions(:,4));
    [M,I] = max(Te_arr);
    time = table2array(arrheniusRates(:,1));
    fig = figure;
    hold on
    for rxn = (str2double(RxnsNum))'
        % plot arrhenius output (1/2) reaction rate
        % arrheniusRates is a table with rates with units of cm^-3 s^-1
        n_arr =  table2array(arrheniusRates(:,rxn+1));
        disp(['Found Arrhenius Rate [1/(cm^3*s^3)] for reaction: ',char(arrheniusRates(:,rxn+1).Properties.VariableNames)])
        % k_arr = n_arr./speciesDens_arr./elecDens_arr;
        plotTitle = char(arrheniusRates(:,rxn+1).Properties.VariableNames);
        % remove the 3 from C3
        plotTitle = erase(plotTitle, '3');
        plot(time, n_arr, 'LineWidth',2, 'DisplayName',plotTitle);
    end
    hold off
    xlim([0 4e-9])
    fig.Children.YAxis.Scale = 'log';
    fontsize(15,'points')
    xlabel('Time [s]')
    ylabel('Reaction Rate [1/cm^3*s]')
    legend
    title('N2(C) Loss Mechanism')
    grid on
    
end