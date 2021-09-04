close all;

% Search in 10x10 grid and create wind farms
props;
cost_params;

g.resetMask();

farm = Windfarm(g, 300, 200, false);
if numel(farm.turbines) > 0
            farm.connect2backbone(g, 50);
            farm.calculate_power();
            farm.calculateCostOnLoc();
            LCOE = farm.LCOE;
end
%%
g.plot_pipes_nav(farm.bbLines);
farm.plot(g);