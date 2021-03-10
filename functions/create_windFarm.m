% Search in 10x10 grid and create wind farms
props;

g.resetMask();

farm = Windfarm(g, 300, 270, false);
farm.connect2backbone(g, 50);

%%
farm.plot(g);
g.plot_pipes();
