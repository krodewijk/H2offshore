% Search in 10x10 grid and create wind farms
close all;

props;
cost_params;

% Load parameters from excel file
[property, costTables] = updatePropsFromExcel("Parameters.xlsm");

% Reset grid mask
g.resetMask();

farm = Windfarm(g, 550, 80, false);
if numel(farm.turbines) > 0
    farm.connect2backbone(g, 50);
    farm.calculate_power();
    farm.calculateCostOnLoc();
    LCOE_onloc = farm.LCOE;
    farm.calculateCost();
    LCOE = farm.LCOE;
end

%% Show results
g.plot_pipes_nav(farm.bbLines);
farm.plot(g);