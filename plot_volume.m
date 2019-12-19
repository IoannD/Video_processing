function [] = plot_volume(VAR)
    % Построение графика выполнения закона сохранения
    % массы от времени
    figure1 = figure('Color',[1 1 1]);
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    ylabel('Volume, %');
    xlabel('Time, s');
    grid(axes1,'on');
    set(axes1,'FontSize',14);
    hold on
    
    referenceVolume = (4/3)*pi*VAR.Radius^3; 
    relativeVolume = 100*(referenceVolume - VAR.volume)/referenceVolume;
    plot(VAR.time,relativeVolume,'.-')
end