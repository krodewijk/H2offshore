global property;

%----------------------
%      SCENARIO          
%----------------------
% What scenario to use?
% 'H2inTurb'     --> Create H2 directly in turbine 
% 'fullElectric' --> Onshore electrolyzer
% 'electricToH2' --> Offshore electrolyzer
property.scenario = 'H2inTurb';


%----------------------
%       GRID          
%----------------------

% Limits of generated map
property.latlim = [50.5, 61.5]; % latitudal (vertical) angle limits in degrees
property.lonlim = [-3.2, 8.2]; % longitudal (horizontal) angle limits in degrees

% Grid resolution
property.resolution = 2.5; % km grid (possible turbine placement every $RESOLUTION kilometers)

% Grid projection system (to convert global geographical latitude and
% longitude coordinates to local km projected coordinate system)
property.projSystem = 23095; % EPSG number (from https://epsg.io/)

% Thresholds to take into account when creating map mask (for possible turbine/objects placements)
property.fishing_threshold = 30; % Minimum average mW fishing hours
property.traffic_threshold = 500;  % Minimum routes/square km per year


%----------------------
%        MAPS          
%----------------------
% Location of Geotif data
property.windHeight = 10; % m (height of measured windspeed)
%property.maps.windSpeed = "GIS_data/windspeed/wind2011_speed.tif";
property.maps.windVariance = "GIS_data/windspeed/wind2011_var.tif";
property.maps.routeDensity = "GIS_data/traffic.tif";
property.maps.bathymetry = "GIS_data/comb_bath.tif";

% Location of shapefiles
property.shapes.natura2000 = "GIS_data/natura2000.shp";
property.shapes.CDDA = "GIS_data/CDDA.shp";
property.shapes.fishing = "GIS_data/Average_MW_Fishing"; % average mW fishing hours
property.shapes.windfarms = "GIS_data/windfarm.shp"; % already existing + under construction windfarms
property.shapes.pipelines = "GIS_data/backbone_simplified/Complete_Backbone_Simplify.shp"; % already existing pipelines (used+unused)
property.shapes.onshoreConnections = "GIS_data/onshoreConnections/Connection Points.shp"; % possible onshore connection points

%----------------------
%      CLUSTERS          
%----------------------

% What connection point to connect to
% 'closest' --> calculate closest connection point from offshore platform
% "NL", "DK", "DE", "NO", "UK", "BE", "FR", "GE" --> Connects to closest connection point in
% corresponding country
%property.shoreConnectionPoint = "NL"; 
property.shoreConnectionPoint = 5;


% clusterMappingOption: How should the cluster layout be?
% "size" --> Cluster is mapped purely by set dimensions (completely fills
% the cluster and puts the hub in centre)
% "power" --> Platform is set at determined coordinates and windfarms are put
% around it until set maximum power is reached. (clusterDim is now used as
% maximum size)

property.clusterMappingOption = "size"; % see explanation above
property.clusterPower = 10; % GW (maximum cluster power)
property.clusterDim = [20, 20]; % [X,Y] in grid coordinates (so X by Y property.resolution meters)
property.bbConnectionPlatform = true;  % Need some kind of transformer/compressor substation near the backbone

%----------------------
%     WINDFARMS          
%----------------------

% Size of windfarm (in grid coordinates if completely full (can be less because of grid mask))
property.farmDim = [10, 10]; % [X, Y] in grid coordinates (so [10,10] is 10 by 10 (=100) turbines at property.resolution km distance)

% Size of windfarm (in total power)
% CAUTION: if > 0 , property.farDim is now used as max size!
property.farmPower = 1; % GW

% Minimum turbines required in a farm
property.farmMinTurbs = 2;

% General baseload of backbone
property.BBbaseload = 0.6; % 60 % baseload


%----------------------
%     TURBINES        
%----------------------

% Turbine parameters
property.turbRating = 12e6; % W (power rating)
property.airDensity = 1.225; % kg/m3 (Air Density)
property.turbRadius = 110; % m (turbine radius)
property.turbMaxWind = 22; % m/s (max rated wind speed)
property.turbMinWind = 4.5; % m/s (min rated windspeed) 
property.turbHubHeight = 260; % m (hub height)
property.turbH2Eff = 0.82; % percent
property.turbWindEff = 0.3; % percentage of wind energy converted to electricity

%----------------------
%     SUBPLATFORMS        
%----------------------
% Platform parameters
property.platform.floating = false; % 

% Electrolyser parameters
property.platform.electrolyserEff = 82; % percent
property.platform.electrolyserOutPress = 30; % Bar

% Transformer parameters
property.platform.transfEff = 95; % percent
property.platform.transOutVolt = 220; % kV

% Compressor parameters
property.platform.comprEff = 80; % percent
property.platform.comprAdiaEff = 90; % percent
property.platform.comprOutPress = 30; %bar

%----------------------
%     CABLES        
%----------------------

property.dielectricLossFactor = 4e-2;

% Turbine --> subStation cable rating
property.turb2subCableVRating = 275; % kV (voltage rating)
property.turb2subCableIRating = 825; % A (current rating)
property.turb2subCableFreq = 50; % Hz
property.turb2subCableA = 1000; % mm2 (cross-section)
property.turb2subCableCap = 0.18; % uF / km (capacitance)
property.turb2subCableRho = 0.0171; % Ohm * mm2/m (copper resistivity)

% substation --> bb cable rating
property.sub2hubCableVRating = 275; % kV (voltage rating)
property.sub2hubCableIRating = 825; % A (current rating)
property.sub2hubCableFreq = 50; % if cable is AC or DC (true = AC, false = DC)
property.sub2hubCableA = 1000; % mm2 (cross-section)
property.sub2hubCableCap = 0.18; % uF / km (capacitance)
property.sub2hubCableRho = 0.0171; % Ohm * mm2/m (copper resistivity)

% Backbone cable rating
property.hub2shoreCableVRating = 525; % kV (voltage rating)
property.hub2shoreCableIRating = 4952; % A (current rating)
property.hub2shoreCableFreq = 0; % if cable is AC or DC (true = AC, false = DC)
property.hub2shoreCableA = 3000; % mm2 (cross-section)
property.hub2shoreCableCap = 0; % uF / km (capacitance)
property.hub2shoreCableRho = 0.0171; % Ohm * mm2/m (copper resistivity)

%----------------------
%      PIPES        
%----------------------
% edit: not sure what values need to be set here... %
property.H2density = 0.08988; % kg/m3 (H2 density @ 288 K, 1bar)
property.H2specEnergy = 39389; % Wh/kg (H2 specific energy)
property.H2compressibility = 1.06; 
property.H2normalGasTemp = 291; % K
property.H2gravity = 0.0696;
property.baseTemp = 288; % K
property.basePress = 1.01325; % bar
property.maxV = 20; % m/s (max velocity)
property.specHeatRatio = 1.41; % Cp / Cv

% Turbine --> hub pipe rating
property.turb2subPipePressure = 30; % bar (pressure rating)
property.turb2subPipeRadius = 0.25; % m (pipe radius)
property.turb2subPipeGasFlow = 339292; % m3/day (gas volumetric flow)
property.turb2subPipeEff = 0.95; % pipeline efficiency
property.turb2subInT = 50; % degC (input temp)
property.turb2suboutT = 10; % degC (output temp)

% substation --> main pipe rating
property.sub2mainPipePressure = 50; % bar (pressure rating)
property.sub2mainPipeRadius = 0.375; % m (pipe radius)
property.sub2mainPipeGasFlow = 763407; % m3/day (gas volumetric flow)
property.sub2mainPipeEff = 0.95; % pipeline efficiency
property.sub2mainInT = 50; % degC (input temp)
property.sub2mainoutT = 10; % degC (output temp)

% main --> shore pipe rating
property.main2bbPipePressure = 80; % bar (pressure rating)
property.main2bbPipeRadius = 0.5; % m (pipe radius)
property.main2bbPipeGasFlow = 1357168; % m3/day (gas volumetric flow)
property.main2bbPipeEff = 0.95; % pipeline efficiency
property.main2bbInT = 50; % degC (input temp)
property.main2bboutT = 10; % degC (output temp)

