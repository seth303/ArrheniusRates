function rateCompareWithEqn(rxnTable,Te,Te_comp,Eqn1,Eqn2)

% fig = figure;
nexttile
    % assemble arrenhius form 
    A = rxnTable.A; Ea = rxnTable.Ea; n = rxnTable.n;

    arrenhiusForm = @(A,n,Ea,Te) A.*Te.^n.*exp(-Ea./Te);
    
    % plottig Wilson
    plot(Te, arrenhiusForm(A,n,Ea,Te), 'LineWidth',2)
    hold on
    
    % plotting comparrison
    plot(Te, Eqn1(Te_comp),'LineWidth',2)

    % plotting comparrison
    plot(Te, Eqn2(Te_comp),'LineWidth',2)
    
    hold off
    % plot options
    % xlabel('Electron Temp [K]')
    % ylabel('Reaction Rate [cm^3/s]')
    grid on
    title(rxnTable.Chemistry)
    % fig.Children.YScale = 'log';
    fig.Children.XAxis.Exponent = 0;
    fontsize(20, "points")
    legend('Wilson' , 'Redondo','Location','best')
    axis padded
end