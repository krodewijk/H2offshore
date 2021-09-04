global costsParams;

% ------------------------- %
%  GLOBAL COST PARAMETERS   %
% ------------------------- %

costsParams.WACC = 0.06;
costsParams.lifetime = 27; % years
costsParams.annuityFactor = 0.0559; 

% -------------------------- %
%  TURBINE COST PARAMETERS   %
% -------------------------- %

costsParams.steelCosts = 500; % EUR/ton
costsParams.embLperDepth = 1.3; 
costsParams.pileA = 6.3; % m2/MW
costsParams.pileDperThickness = 49; % m
costsParams.thickIncr = 0.07; % %/5m depth
costsParams.steelWeight = 7.85; % ton/m3
costsParams.transitionLength = 26; % meter

costsParams.foundInstallation = 0.3; % M-EUR/1000 ton
costsParams.foundMinSteel = 0.699; % M-EUR/1000 ton

costsParams.turbCosts = 0.65; % M-EUR/MW
costsParams.turbCostsExtra = 0.4; % M-EUR/MW
costsParams.turbInstallation = 0.041; % M-EUR/MW
costsParams.turbWeight = 140; % ton/MW


% -------------------------- %
%    CABLE COST PARAMETERS   %
% -------------------------- %

costsParams.turb2subCable = 1; % M-EUR/km
costsParams.hub2shoreCable = 3.1; % M-EUR/km

% -------------------------- %
%    PIPE COST PARAMETERS   %
% -------------------------- %

costsParams.turb2subPipe = 0.1; % M-EUR/km
costsParams.hub2shorePipe = 1; % M-EUR/km


% --------------------------- %
%  PLATFORM COST PARAMETERS   %
% --------------------------- %

costsParams.platform = 400; % M-EUR
costsParams.H2ProductionCAPEX = 150; % EUR/kW
costsParams.comprCost = 2; % M-EUR/MW
costsParams.transfCost = 176.75; % M-EUR/GW

% --------------------------- %
%     OTHER COST PARAMETERS   %
% --------------------------- %

costsParams.management = 119; % M-EUR / GW
costsParams.insurance = 0.10; % of total CAPEX
costsParams.maintenance = 0.02; % of CAPEX (equip + instal)
