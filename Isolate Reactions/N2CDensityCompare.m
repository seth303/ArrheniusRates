clc; clear; close all;
Delays0s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 0s.dat');
Delays2s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 2 ns.dat');
Delays4s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 4 ns.dat');
Delays6s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 6 ns.dat');

% ID import options, and set them
arrheniusRatesBranch = gitrepo('Ivanov/');
speciesListFilepath =      strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_species_list.txt');
arrDensitiesFilePath =     strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_300K_FIELD0.2/qt_densities.txt');
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

%% plot the time resolved N2(C)
figure;
plot(Delays0s.("Wavelength [nm]"),Delays0s.("Intensity [au]"),'LineWidth',2)
hold on
plot(Delays2s.("Wavelength [nm]"),Delays2s.("Intensity [au]"),'LineWidth',2)
plot(Delays4s.("Wavelength [nm]"),Delays4s.("Intensity [au]"),'LineWidth',2)
plot(Delays6s.("Wavelength [nm]"),Delays6s.("Intensity [au]"),'LineWidth',2)
hold off
legend('delay = 0ns','delay = 2ns','delay = 4ns','delay = 6ns')
xlabel('Wavelength [nm]')
ylabel('Intensity [a.u.]')
title('Emission Spectra of N2(C)')
fontsize(15, 'points')
grid on

%% compute the relative population at each level
% close all
n_N2C_Ons = sum(Delays0s.("Intensity [au]"))/max(Delays0s.("Wavelength [nm]"));
n_N2C_2ns = sum(Delays2s.("Intensity [au]"))/max(Delays2s.("Wavelength [nm]"));
n_N2C_4ns = sum(Delays4s.("Intensity [au]"))/max(Delays4s.("Wavelength [nm]"));
n_N2C_6ns = sum(Delays6s.("Intensity [au]"))/max(Delays6s.("Wavelength [nm]"));

n_N2C = [n_N2C_Ons n_N2C_2ns n_N2C_4ns n_N2C_6ns];
n_N2C = rescale(n_N2C);
time_Experimental = (0:2:6)*1e-9;
fig = figure;
plot(time_Experimental,n_N2C, '-*','LineWidth',2)
hold on
% plotting experimental data
timeArrhenius = table2array(arrhenius_Densities(:,1));
Individual = table2array(arrhenius_Densities(:,14));
% normalize Individual
Individual = rescale(Individual);
plot(timeArrhenius,Individual,'-*','LineWidth',2)

% plot options
hold off
fig.Children.XAxis.Exponent = -9;
legend('Experimental','Model')
xlabel('Time [s]')
ylabel('Relative Population [a.u.]')
title('Emission Spectra of N2(C)')
fontsize(15, 'points')
grid on
xlim([0 10e-9])

