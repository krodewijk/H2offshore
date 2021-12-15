folder = 'results/08-Dec-2021/';

% Load all files
fullElectric = load([folder, 'fullElectric']);
%H2inTurb = load([folder, 'H2inTurb']);

% filter out all zeros and values higher than 300
fullElectric.LCOE(fullElectric.LCOE == 0 | fullElectric.LCOE > 300) = -1;
fullElectric.LCOE_onloc(fullElectric.LCOE_onloc == 0 | fullElectric.LCOE_onloc > 300) = -1;

% H2inTurb.LCOE(H2inTurb.LCOE == 0 | H2inTurb.LCOE > 300) = -1;
% H2inTurb.LCOE_onloc(H2inTurb.LCOE_onloc == 0 | H2inTurb.LCOE_onloc > 300) = -1;
% 
% H2inTurb.LCOH(H2inTurb.LCOH == 0 | H2inTurb.LCOH > 10) = -1;
% H2inTurb.LCOH_onloc(H2inTurb.LCOH_onloc == 0 | H2inTurb.LCOH_onloc > 10) = -1;
% 
% H2inTurb.floating(H2inTurb.LCOE == 0 | H2inTurb.LCOE > 300) = -1;

%% Save data to tiff files
geotiffwrite([folder, 'LCOE_fullElectric'], fullElectric.LCOE, g.cellReference);
geotiffwrite([folder, 'LCOE_fullElectricOnLoc'], fullElectric.LCOE_onloc, g.cellReference);

%%
geotiffwrite([folder, 'LCOE_H2inTurb'], H2inTurb.LCOE, g.cellReference);
geotiffwrite([folder, 'LCOE_H2inTurbOnLoc'], H2inTurb.LCOE_onloc, g.cellReference);

geotiffwrite([folder, 'LCOH_H2inTurb'], H2inTurb.LCOH, g.cellReference);
geotiffwrite([folder, 'LCOH_H2inTurbOnLoc'], H2inTurb.LCOH_onloc, g.cellReference);

geotiffwrite([folder, 'floating'], floating, g.cellReference);
