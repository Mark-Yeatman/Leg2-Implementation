%Plots several figures/subplots in one tabbed window

set(0,'DefaultAxesXGrid','on');
set(0,'DefaultAxesYGrid','on');
set(0,'DefaultFigureWindowStyle','docked');

Analysis;

plotTorques;
plotPosVel;
plotPower;
plotKPBC_Data;
plotLoadCell;

set(0,'DefaultAxesXGrid','off');
set(0,'DefaultAxesYGrid','off');
set(0,'DefaultFigureWindowStyle','default');