% clc; clear; 
clc; close all;
load ../MASTER_BOLSIG_Rates_Data.mat
n = 67;


reactionName = sprintf('Reaction_%d',n);
epsilon = data.(reactionName).data(:,1);
sigma = data.(reactionName).data(:,2);

% remove zeros from sigma
% r = sigma>0;
% sigma = sigma(r);
% epsilon = epsilon(r);
% % normalize sigma
% sigma = sigma/max(sigma);
% sigma = sigma*1e22;



% 
% [paramEsts,paramCIs] = gevfit(sigma);
% 
% kMLE = paramEsts(1)        % Shape parameter
% sigmaMLE = paramEsts(2)    % Scale parameter
% muMLE = paramEsts(3)       % Location parameter
% 
% kCI = paramCIs(:,1)
% sigmaCI = paramCIs(:,2)
% muCI = paramCIs(:,3)
% 
% [nll,acov] = gevlike(paramEsts,sigma);
% paramSEs = sqrt(diag(acov))
% 
% lowerBnd = muMLE-sigmaMLE./kMLE;
% 
% ymax = 1.1*max(sigma);
% bins = floor(lowerBnd):ceil(ymax);
% h = bar(bins,histc(sigma,bins)/length(sigma),'histc');
% h.FaceColor = [.9 .9 .9];
% ygrid = linspace(lowerBnd,ymax,100);
% line(ygrid,gevpdf(ygrid,kMLE,sigmaMLE,muMLE));
% xlabel('Block Maximum');
% ylabel('Probability Density');
% xlim([lowerBnd ymax]);
% 
% 


plot(epsilon,sigma,'o');
xlabel('Energy [eV]');
ylabel('Cross Section [m^2]');

pd =fitdist((epsilon), 'InverseGaussian')

hold on
x = linspace(0,max(epsilon));
plot(x,pdf(pd,x)*max(sigma)*mean(epsilon),'LineWidth',2)
% xlim([0 150])
% ylim([0 2e-22])
hold off
title(data.(reactionName).chemistry)

%%

modelFun =  @(p,x) p(3) .* (x./p(1)).^(p(2)-1) .* exp(-(x./p(1)).^p(2));
startingVals = [18 5 20];
nlModel = fitnlm(epsilon,sigma,modelFun,startingVals);

xgrid = linspace(0,200,200)';
line(xgrid,predict(nlModel,xgrid),'Color','r');
% 
% nlModel2 = fitnlm(epsilon,log(sigma),@(p,x) log(modelFun(p,x)),startingVals);
% 
% line(xgrid,exp(predict(nlModel2,xgrid)),'Color',[0 .5 0],'LineStyle','--');
% legend({'Raw Data','Additive Errors Model','Multiplicative Errors Model'});
