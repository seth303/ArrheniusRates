function plotLogLog(x,y,plotTitle,XLabel,YLabel)
    % Plot_LogLog = figure;
    loglog(x,y, '-*')
    grid on
    title(plotTitle);
    xlabel(XLabel);
    ylabel(YLabel);
    fontsize(15, "points")
end