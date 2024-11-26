clc; clear; %close all;
load ../MASTER_BOLSIG_Rates_Data.mat
load ../MASTER_reactionTable.mat


% find Bolsig RXN's
[C, idxA,idxB,bolsigRxnsName] =  findBolsigRxns(reactionTable);
% sort
[idxB,sortIdx] = sort(idxB,'ascend');
idxA = idxA(sortIdx);
C = C(sortIdx);
fig = figure;
fig.WindowState = 'maximized';
tL = tiledlayout;
%% Constants
Te = 3000:100:60000;
% Boltzmann constant in eV/K
kB_eV = 8.617333262145e-5;
Te_eV = Te*kB_eV;

%% I-1 in Redondo -> 26th BOLSIG in Wilson
rxn =reactionTable(idxA(26),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,1.03e-8,0.5,14.3);

%% I-3 in Redondo -> 28th BOLSIG in Wilson
rxn =reactionTable(idxA(28),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,7.07e-9,0.5,13.1);

%% I-12 in Redondo -> 29th BOLSIG in Wilson
rxn =reactionTable(idxA(29),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,9.0e-10,2,11.6);

%% A1 in Redondo -> 31st BOLSIG in Wilson
rxn =reactionTable(idxA(31),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,1.07e-9,-1.391,6.26);

%% X1 in Redondo -> 13th BOLSIG in Wilson
rxn =reactionTable(idxA(13),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,1.7e-9,0,3.1);

%% X2 in Redondo -> 21st BOLSIG in Wilson
rxn =reactionTable(idxA(21),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,4.5e-9,0,2.29);

%% O3->O+O2^- in Schneider -> 34th BOLSIG in Wilson
rxn =reactionTable(idxA(34),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,1.1e-15*(100^3),0,0); % Schneider reports rates in terms of m3/s
legend('Wilson','Schneider')

%% O3->O^- + O2 in Schneider -> 33rd BOLSIG in Wilson
rxn =reactionTable(idxA(33),:);
disp(rxn)
rateCompare(rxn,Te,Te_eV,1.1e-17*(100^3),0,0); % Schneider reports rates in terms of m3/s
legend('Wilson','Schneider')

%% N2->N2^+ in Ali & Schneider -> 27th BOLSIG in Wilson
rxn =reactionTable(idxA(27),:);
disp(rxn)
N2_Ion_Ali = @(Te) 1E-9./sqrt(Te) .* (-4.98 + 9.67.*Te + 0.49*Te.^2) .* exp(-15.58./Te);
% P in torr, Tg in K
Tg = 300; %K
P = 760; % torr
N2_Ion_Sch = @(Te) (P*(300/Tg)*(1E7./sqrt(Te))).*((14.093 + 27.366.*Te + 1.386.*Te.^2).*(exp(-15.58./Te)));

rateCompareWithEqn(rxn,Te,Te_eV,N2_Ion_Ali,N2_Ion_Sch);
legend('Wilson','Ali', 'Schneider')

%% O2->O2^+ in Ali & Schneider -> 28th BOLSIG in Wilson
rxn =reactionTable(idxA(28),:);
disp(rxn)
O2_Ion_Ali = @(Te) 1E-9./sqrt(Te) .* (-2.02 + 3.39.*Te + 0.64*Te.^2) .* exp(-12.06./Te);
O2_Ion_Sch = @(Te) (P*(300/Tg)*(1E7./sqrt(Te))).*((1.43 + 2.47.*Te + 0.456.*Te.^2).*(exp(-12.06./Te)));


rateCompareWithEqn(rxn,Te,Te_eV,O2_Ion_Ali,O2_Ion_Sch);
legend('Wilson','Ali', 'Schneider')

% plot options
xlabel(tL, 'Electron Temp [K]', 'FontSize',25)
ylabel(tL, 'Reaction Rate [cm^3/s]', 'FontSize',25)
for i = 1:7
    if isgraphics(tL.Children(i),'Axes')
    tL.Children(i).YScale = 'log';
    % tL.Children(i).XAxis.Exponent = 0;
    end
end