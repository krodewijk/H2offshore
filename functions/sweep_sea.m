close all;

global property;

% Search in entire north sea and create wind farms
props;
cost_params;

g.resetMask();

xIntrin_size = size(g.lat, 1);
yIntrin_size = size(g.lat, 2);
x_step = property.farmDim(1);
y_step = property.farmDim(2);

LCOE = zeros(size(g.lat));
LCOH = zeros(size(g.lat));

fprintf(1,' Progress:    ');
for x = 1:x_step:xIntrin_size-y_step
    for y = 1:y_step:yIntrin_size-y_step
        farm = Windfarm(g,x,y,false);
        if numel(farm.turbines) > 0
            farm.connect2backbone(g, 50);
            farm.calculate_power();
            %farm.calculateCostOnLoc();
            farm.calculateCost();
            LCOE(x:(x+x_step), y:(y+y_step)) = farm.LCOE;
            LCOH(x:(x+x_step), y:(y+y_step)) = farm.LCOH;
        end
        %display(num2str(y))
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(x/xIntrin_size)*100);
end
fprintf(1,'\n');

save("results3_H2inTurb", "LCOE", "LCOH");
