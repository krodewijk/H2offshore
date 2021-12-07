close all;

props;
cost_params;
global costTables;

fileName = "Parameters.xlsm";
% dimensions of North Sea are approx 470 x 590
x = 140;
y = 370;

% Load parameters from excel file
[property, costTables] = updatePropsFromExcel(fileName);
property.scenario = "fullElectric";

% Reset grid mask
g.resetMask();

farm = Windfarm(g, y, x, false);
if numel(farm.turbines) > 0
    farm.connect2backbone(g, 50);
    farm.calculate_power();
    farm.calculateCost(g);
end

%% Write results to excel sheet
writeResultsToExcel(farm);

%% Do same but now for hydrogen scenario
property.scenario = "H2inTurb";
% Reset grid mask
g.resetMask();

farm = Windfarm(g, y, x, false);
if numel(farm.turbines) > 0
    farm.connect2backbone(g, 50);
    farm.calculate_power();
    farm.calculateCost(g);
end

%% Write results to excel sheet
writeResultsToExcel(farm);

%% Show results
g.plot_pipes_nav(farm.bbLines);
farm.plot(g);