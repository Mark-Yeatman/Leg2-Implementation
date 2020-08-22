function animate(drawfunc, x, t, links)
% Written by Mark Yeatman

    % Initialize the movie
    FramesPerSec = 30;

    mov = VideoWriter('Test');
    mov.FrameRate = FramesPerSec;
    mov.Quality = 100;
    open(mov)

    % Create figure
    screenSize = get(0, 'ScreenSize');
    width = screenSize(3)*1/2;
    left = screenSize(3)/4;
    bottom = screenSize(4)/4;
    height = screenSize(4)*1/2;
    fig = figure('Position', [left bottom width height],...
                 'Color','w',...
                 'DoubleBuffer','on');
    set(gca,'NextPlot','replace','Visible','off')
   
    for i=1:length(t)                                     
        drawfunc(x(i,:)',links);
        pause(0.01)       
        
         % Add movie frames?
        F = getframe(gcf);
        writeVideo(mov,F);
    
    end
    

end