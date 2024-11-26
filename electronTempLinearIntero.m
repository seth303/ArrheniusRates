clc; clear; close all;
load('MASTER_electronTempTimeEvolution.mat');
xq = [0,logspace(-12,-7,200)];
vq1 = interp1(time,T_Te,xq);
plot(time,T_Te,'o', 'MarkerSize',10)
hold on


plot(xq,vq1,'^', 'MarkerSize',4);
legend('Original','Linear Interp')
fontsize(15,'points')
grid on
xlabel('Time [s]')
ylabel('Electron Temp [K]')
title('Logspace Electron Temperature')