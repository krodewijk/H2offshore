close all;

global property;
global costTables;

fileName = "Parameters.xlsm";

% Search in entire north sea and create wind farms
props;
cost_params;

% Load parameters from excel file
[property, costTables] = updatePropsFromExcel(fileName);
property.scenario = "fullElectric";

g.resetMask();

xIntrin_size = size(g.lat, 1);
yIntrin_size = size(g.lat, 2);
x_step = property.farmDim(1);
y_step = property.farmDim(2);

LCOE = zeros(size(g.lat));
LCOE_onloc = zeros(size(g.lat));
floating = zeros(size(g.lat));

fprintf(1,' Progress:    ');
for x = 1:x_step:xIntrin_size-y_step
    for y = 1:y_step:yIntrin_size-y_step
        farm = Windfarm(g,x,y,false);
        if numel(farm.turbines) > 0
            farm.connect2backbone(g, 50);
            farm.calculate_power();
            farm.calculateCost(g);
            LCOE(x:(x+x_step), y:(y+y_step)) = farm.LCOEOnshore;
            LCOE_onloc(x:(x+x_step), y:(y+y_step)) = farm.LCOEOnLoc;
            floating(x:(x+x_step), y:(y+y_step)) = farm.turbines(1).floating;
        end
        %display(num2str(y))
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(x/xIntrin_size)*100);
end
fprintf(1,'\n');

date = datestr(datetime('now'),1);
mkdir(['results/', date]);
save(['results/', date, '/', farm.scenario], "LCOE", "LCOE_onloc");
