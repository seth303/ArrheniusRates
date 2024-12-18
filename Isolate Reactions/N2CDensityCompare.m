clc; clear; %close all;
Delays0s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 0s.dat');
Delays2s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 2 ns.dat');
Delays4s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 4 ns.dat');
Delays6s = importSpectra('380nm/fs-Filament Spectra, Exposure 2 ns, Delay 6 ns.dat');

% ID import options, and set them
arrheniusRatesBranch = gitrepo('Ivanov/');
speciesListFilepath =      strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_species_list.txt');
arrDensitiesFilePath =     strcat(arrheniusRatesBranch.WorkingFolder,'/Results/NRP_300um_500K_FIELD0.2/qt_densities.txt');
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
close all
n_N2C_Ons = sum(Delays0s.("Intensity [au]"))/max(Delays0s.("Wavelength [nm]"));
n_N2C_2ns = sum(Delays2s.("Intensity [au]"))/max(Delays2s.("Wavelength [nm]"));
n_N2C_4ns = sum(Delays4s.("Intensity [au]"))/max(Delays4s.("Wavelength [nm]"));
n_N2C_6ns = sum(Delays6s.("Intensity [au]"))/max(Delays6s.("Wavelength [nm]"));

n_N2C = [n_N2C_Ons n_N2C_2ns n_N2C_4ns n_N2C_6ns];
n_N2C = rescale(n_N2C);
time_Experimental = (0:2:6)*1e-9;
figure1 = figure;
plot(time_Experimental,n_N2C, '--*','MarkerFaceColor',[0 0.447058823529412 0.741176470588235],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',9,...
    'Marker','o',...
    'LineWidth',3,...
    'LineStyle','--');
hold on
% plotting experimental data
timeArrhenius = table2array(arrhenius_Densities(:,1));
Individual = table2array(arrhenius_Densities(:,14));
% normalize Individual
Individual = rescale(Individual);
plot(timeArrhenius,Individual,'-','LineWidth',4)

% plot options
hold off
figure1.Children.XAxis.Exponent = -9;
legend('Experimental','Model')
xlabel('Time [s]')
ylabel('Relative Population [a.u.]')
title('Emission Spectra of N2(C)')
fontsize(15, 'points')
grid on
xlim([0 10e-9])

%%
% to find the first data point, I need to integrate the model curve until
% it reaches a value equal to the first experimental data point, and then
% go backwards 2ns-X to find the place to put the model curve. Then, I need
% to simulate gate widths of 2ns to show that when the model is integrated,
% the values approach that of the experiment
gateWidth = 2e-9;
% cumulative sum of model
Q = cumsum(Individual);
Q = rescale(Q);
% first datapoint
experiment = n_N2C(1);
% find where they intersect
[~,closestIndex_frame1] = min(abs(experiment-Q));
% time of first datapoint
realDelay0 = timeArrhenius(closestIndex_frame1);

% find time of 2ns after real time 0
realGateWidth2 = realDelay0  + gateWidth;
% find where they intersect
[~,closestIndex_frame2] = min(abs(timeArrhenius-realGateWidth2));
% time of first datapoint
realDelay2 = timeArrhenius(closestIndex_frame2);

% find time of 2ns after real time 0
realGateWidth4 = realDelay2  + gateWidth;
% find where they intersect
[~,closestIndex_frame3] = min(abs(timeArrhenius-realGateWidth4));
% time of first datapoint
realDelay4 = timeArrhenius(closestIndex_frame3);

% find time of 2ns after real time 0
realGateWidth6 = realDelay4  + gateWidth;
% find where they intersect
[~,closestIndex_frame4] = min(abs(timeArrhenius-realGateWidth6));
% time of first datapoint
realDelay6 = timeArrhenius(closestIndex_frame4);

rectplot = @(x1,x2) rectangle('Position',[x1 0 (x2-x1) 1], 'FaceColor','#7E2F8E','FaceAlpha',0.3,'Curvature',[.1 .1]);


hold on
frame1 = rectplot(realDelay0-gateWidth, realDelay0);
frame2 = rectplot(realDelay2-gateWidth, realDelay2);
frame3 = rectplot(realDelay4-gateWidth, realDelay4);
frame4 = rectplot(realDelay6-gateWidth, realDelay6);
xlim([-2.1e-9 7e-9])
hold off

frame1Integral = trapz(timeArrhenius(1:closestIndex_frame1),Individual(1:closestIndex_frame1));
frame2Integral = trapz(timeArrhenius(closestIndex_frame1:closestIndex_frame2),Individual(closestIndex_frame1:closestIndex_frame2));
frame3Integral = trapz(timeArrhenius(closestIndex_frame2:closestIndex_frame3),Individual(closestIndex_frame2:closestIndex_frame3));
frame4Integral = trapz(timeArrhenius(closestIndex_frame3:closestIndex_frame4),Individual(closestIndex_frame3:closestIndex_frame4));

% frame1Integral = sum(Individual(1:closestIndex_frame1));
% frame2Integral = sum(Individual(closestIndex_frame1:closestIndex_frame2));
% frame3Integral = sum(Individual(closestIndex_frame2:closestIndex_frame3));
% frame4Integral = sum(Individual(closestIndex_frame3:closestIndex_frame4));

framesIntegral = [frame1Integral frame2Integral frame3Integral frame4Integral];
framesIntegral = rescale(framesIntegral);
hold on
plot(time_Experimental,framesIntegral,    'MarkerFaceColor',[0.929411764705882 0.694117647058824 0.125490196078431],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',9,...
    'Marker','o',...
    'LineWidth',3,...
    'LineStyle','--',...
    'DisplayName','Model - Exposure Time Corrected');
% fill([0 0 0 0], [0 0 0 0],[0.4940 0.1840 0.5560], 'FaceAlpha',0.3, 'DisplayName','Frames')

annotation(figure1,'textbox',[0.2 0.5 0.08 0.08],'Color',[1 1 1],...
    'String','Frame 1',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'FaceAlpha',0.8,...
    'BackgroundColor',[0.494117647058824 0.184313725490196 0.556862745098039]);

annotation(figure1,'textbox',[0.38 0.5 0.08 0.08],'Color',[1 1 1],...
    'String','Frame 2',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'FaceAlpha',0.8,...
    'BackgroundColor',[0.494117647058824 0.184313725490196 0.556862745098039]);

annotation(figure1,'textbox',[0.55 0.5 0.08 0.08],'Color',[1 1 1],...
    'String','Frame 3',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'FaceAlpha',0.8,...
    'BackgroundColor',[0.494117647058824 0.184313725490196 0.556862745098039]);

annotation(figure1,'textbox',[0.7 0.5 0.08 0.08],'Color',[1 1 1],...
    'String','Frame 4',...
    'FontSize',24,...
    'FitBoxToText','off',...
    'FaceAlpha',0.8,...
    'BackgroundColor',[0.494117647058824 0.184313725490196 0.556862745098039]);

