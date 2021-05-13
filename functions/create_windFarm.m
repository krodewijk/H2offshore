% Search in 10x10 grid and create wind farms
props;
cost_params;

g.resetMask();

farm = Windfarm(g, 490, 370, false);
farm.connect2backbone(g, 50);
farm.calculate_power();
farm.calculateCost();
%%
g.plot_pipes_nav(farm.bbLines);
farm.plot(g);