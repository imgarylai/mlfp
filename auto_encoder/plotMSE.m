function plotMSE(prefix, type, mse)
    f_name = sprintf('%s-%s', prefix, type);
    [row, column] = size(mse);
    meanMse = mean(mse);
    stdMse = std(mse);
    fig = figure;
    hold on;
    title(f_name);
    xlabel('Numbers of neuron in the hidden layer.');
    ylabel('MSE');
    H1 = plot(1:column, meanMse, '-', 'lineWidth', 2);
    H2 = plot(1:column, [meanMse - stdMse; meanMse + stdMse],':', 'lineWidth', 2);
    legend([H1,H2(1)],'mse','mse std', 'Location', 'Northeast');
    hold off;
    saveas(fig, f_name, 'fig');
    saveas(fig, f_name, 'png');
end