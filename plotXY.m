function plotXY(x,y, plotTitle, XLabel, YLabel)
    plot(x, y, '-*')
    grid on
    title(plotTitle);
    xlabel(XLabel);
    ylabel(YLabel);
    fontsize(15, "points")
end