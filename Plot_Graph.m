function Plot_Graph(x,y,plotTitle,xLabel,yLabel)
    figure
    plot(x, y, '-m');
    title(plotTitle);
    xlabel(xLabel);
    ylabel(yLabel);
    grid on;
end