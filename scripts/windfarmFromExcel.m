% Search in 10x10 grid and create wind farms
close all;

props;
cost_params;
global costTables;

fileName = "Parameters.xlsm";

% Load parameters from excel file
[property, costTables] = updatePropsFromExcel(fileName);

% Reset grid mask
g.resetMask();

farm = Windfarm(g, 550, 100, false);
if numel(farm.turbines) > 0
    farm.connect2backbone(g, 50);
    farm.calculate_power();
    farm.calculateCost();
end

%% Write results to excel sheet
writeResultsToExcel(farm);

%% Show results
g.plot_pipes_nav(farm.bbLines);
farm.plot(g);