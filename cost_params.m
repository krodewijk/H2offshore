global costsParams;

% ------------------------- %
%  GLOBAL COST PARAMETERS   %
% ------------------------- %

costsParams.WACC = 0.032;
costsParams.lifetime = 27; % years
costsParams.annuityFactor = 0.0559; 

% -------------------------- %
%  TURBINE COST PARAMETERS   %
% -------------------------- %

costsParams.steelCosts = 500; % EUR/ton
costsParams.embLperDepth = 1.3; 
costsParams.pileA = 6.3; % m2/MW
costsParams.pileDperThickness = 70; % m
costsParams.pileDperThickness_2050reduction = 0.3;
costsParams.thickIncr = 0.07; % %/5m depth
costsParams.steelWeight = 7.85; % ton/m3
costsParams.transitionLength = 26; % meter

costsParams.foundInstallation = 0.428; % M-EUR/1000 ton
costsParams.foundInstallation_2050reduction = 0.3;
costsParams.foundMinSteel = 0.699; % M-EUR/1000 ton

costsParams.turbCosts = 0.760; % M-EUR/MW
costsParams.turbCosts_2050reduction = 0.15;
costsParams.turbCostsExtra = 0.4; % M-EUR/MW
costsParams.turbInstallation = 0.0580; % M-EUR/MW
costsParams.turbInstallation_2050reduction = 0.3;
costsParams.turbWeight = 140; % ton/MW


% -------------------------- %
%    CABLE COST PARAMETERS   %
% -------------------------- %

costsParams.cableInstallation = 1.37; % M-EUR/km
costsParams.turb2subCosts = 0.33; % M-EUR/km
costsParams.hub2shoreCosts = 2.53; % M-EUR/km

% -------------------------- %
%    PIPE COST PARAMETERS   %
% -------------------------- %

costsParams.turb2subPipe = 1; % M-EUR/km
costsParams.hub2shorePipe = 2; % M-EUR/km


% --------------------------- %
%  PLATFORM COST PARAMETERS   %
% --------------------------- %

costsParams.platformEquipm = 145; % M-EUR
costsParams.platformInstallation = 41; % M-EUR
costsParams.H2ProductionCAPEX = 250; % EUR/kW

% --------------------------- %
%     OTHER COST PARAMETERS   %
% --------------------------- %

costsParams.management = 140; % M-EUR / GW
costsParams.management_2050reduction = 0.15;
costsParams.insurance = 0.10; % of total CAPEX
costsParams.maintenance = 0.02; % of CAPEX (equip + instal)
