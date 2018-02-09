function [fig] = plotS2(f,varargin)
% plotS2 - plot a signal on the sphere
    numSignals = size(f,3);
    fig = gcf;
    X = cell(1,3);
    [X{:}] = sphere(40);
    lightGrey = 0.8*[1 1  1]; % It looks better if the lines are lighter
    surface(X{:},'FaceColor', 'none','EdgeColor',lightGrey,'LineWidth',.1)
    hold on
    for i=1:numSignals
        plot3(f(1,:),f(2,:),f(3,:),varargin{:}) 
    end
    hold off
    axis off;
end

