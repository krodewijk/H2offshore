close all;

props;
cost_params;
global costTables;

[property, costTables] = updatePropsFromExcel("Parameters.xlsm");

g.resetMask();

farm = Windfarm(g, 550, 80, false);
if numel(farm.turbines) > 0
    farm.connect2backbone(g, 50);
    farm.calculate_power();
    farm.calculateCost();
    LCOE_onshore = farm.LCOEOnshore;
    LCOE_onLoc = farm.LCOEOnLoc;
end
%%
g.plot_pipes_nav(farm.bbLines);
farm.plot(g);