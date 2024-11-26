%% Import Table
clc; clear; %close all;
% define the git repo
mainBranch = gitrepo('main/');
% arrheniusRatesBranch = gitrepo('arrheniusRates/');
% arrheniusRatesBranch = gitrepo('fsPostFilament/');
arrheniusRatesBranch = gitrepo('Ivanov/');


pull(mainBranch)
pull(arrheniusRatesBranch)
arrheniusRatesBranch.LastCommit.Message


currBranch = arrheniusRatesBranch.CurrentBranch.Name;
% define filepaths
if contains(arrheniusRatesBranch.WorkingFolder,'arrheniusRates')
    arrDensitiesFilePath =     strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt');
    speciesListFilepath =      strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_species_list.txt');
    bolsigDensitiesFilePath =  strcat(mainBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt');
    timescale = 'ns';
end
% fs
if contains(arrheniusRatesBranch.WorkingFolder,'fsPostFilament') || contains(arrheniusRatesBranch.WorkingFolder,'Ivanov') 
    arrDensitiesFilePath =     strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_densities.txt');
    speciesListFilepath =      strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_species_list.txt');
    bolsigDensitiesFilePath =  strcat(mainBranch.WorkingFolder,'/Results/NRP_300um_1200K_FIELD0.02/qt_densities.txt');
    timescale = 'fs';
end

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