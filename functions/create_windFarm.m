% Search in 10x10 grid and create wind farms
props;

g.resetMask();

farm = Windfarm(g, 490, 370, true);
farm.connect2backbone(g, 50);
farm.calculate_power();
%%
farm.plot(g);
g.plot_pipes();
