function plotXYY(xleft,yleft, xright, yright, PlotTitle, XLabel, YLeftLabel, YRightLabel)
    yyaxis left
    plot(xleft,yleft, '*-')
    hold on
    yyaxis right
    plot(xright,yright, '-*', 'Color',"#D95319")
    grid on
    hold off

    fontsize(15, "points")

    title(PlotTitle)
    xlabel(XLabel)
    ylabel(YLeftLabel)
    yyaxis right
    ylabel(YRightLabel)

end