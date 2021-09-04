classdef Turbine < handle
    % Class decribing a single turbine
    
    properties
        rating; % W (power rating)
        density = 1.225; % kg/m3
        A; % m^2
        max_wind; % m/s
        min_wind; % m/s
        windEff;
        
        depth; % m
        wind; % m/s
        windVar;
        windDist;
        windDistBins;
        power; % W (actually generated power)
        hubHeight; % m
        H2out; % true for inturbine electrolyser, false for electric output
        H2Eff;
        
        xIntrin; % with reference to defined grid 
        yIntrin;
        
        foundationWeight;
        OPEX;
        CAPEX; % M-EUR
        
        node; % node in grid graph
    end
    
    methods
        function obj = Turbine(grid, xIntrin, yIntrin)
            global property;
            obj.rating = property.turbRating;
            obj.density = property.airDensity; % kg/m3
            obj.A = pi * property.turbRadius^2;
            obj.max_wind = property.turbMaxWind; % m/s
            obj.hubHeight = property.turbHubHeight;
            obj.H2Eff = property.turbH2Eff;
            obj.windEff = property.turbWindEff;
            obj.min_wind = property.turbMinWind; % m/s
            
            if property.scenario == "H2inTurb"
                obj.H2out = true;
            else
                obj.H2out = false;
            end
            
            % Constructor
            obj.xIntrin = xIntrin;
            obj.yIntrin = yIntrin;
            
            obj.depth = grid.data.bath(xIntrin, yIntrin);
            obj.wind = grid.data.wind(xIntrin, yIntrin);
            obj.windVar = grid.data.windVar(xIntrin, yIntrin);
            obj.windDist = grid.data.windDist(xIntrin, yIntrin, :);
            obj.windDistBins = grid.data.windDistBins;
            obj.node = grid.intrin2node(xIntrin, yIntrin);
            
            obj.calculateExpectedPower();
        end
        
        function windHub = translate_wind(obj, wind)
            global property;
            wind60 = wind .* log(60/0.0001) ./ log(property.windHeight/0.0001);
            windHub = wind60 .* (obj.hubHeight ./ 60) .^ 0.115;
        end
        
        function calculateExpectedPower(obj)
            x = obj.windDistBins;
            
            % Translate wind bins
            x = obj.translate_wind(x);
            
            %fx = 1 / (std*sqrt(pi * 2)) .* exp(-0.5 .* ((x-obj.wind) ./ std).^2); % Wind distribution
            
            % Import histogram with normalized frequencies
            hist = squeeze(obj.windDist)';
            
            % Make from histogram a probability density function
            dist = (1/trapz(x,hist)) .* hist;
            
            % Create PV curve
            curve = zeros(1, numel(x));
            for i = 1:numel(x)
                curve(i) = obj.calculateWindPower(x(i));
            end
            
            % Multiply random distribution with curve and integrate
            expected = trapz(x, dist.*curve);
            
            if obj.H2out
                obj.power = obj.H2Eff * expected;
            else
                obj.power = expected;
            end
        end
        
        function plot_pvCurve(obj)
            x = 0.1:0.1:25;
            y = zeros(1, numel(x));
            for i = 1:numel(x)
                y(i) = obj.calculateWindPower(x(i));
            end
            plot(x, y);
            title("PV-Curve Turbine");
            xlabel("Wind (m/s)");
            ylabel("Power (W)");
            ylim([0, 30e6]);
        end
        
        function turb_power = calculateWindPower(obj, wind)
            if (wind < obj.max_wind) && (wind > obj.min_wind)
                % Calculate wind power
                wind_power = 0.5 * obj.density * obj.A * wind^3;

                % Assume efficiency of ?? %
                turb_power = obj.windEff * wind_power;

                % Adjust for rating
                turb_power(turb_power > obj.rating) = obj.rating;
            else
                turb_power = 0;
            end
        end
        
        function cost = calculateCost(obj)
            global costsParams;
            
            % Equipment costs (in M-EUR)
            
            turbEq = (obj.rating / 1e6) * costsParams.turbCosts;
            
%             if obj.H2out
%                 turbEq = turbEq + (obj.rating / 1000) * (costsParams.H2ProductionCAPEX / 1e6);
%             end
            
            warranty = 0.53 * turbEq;
            
            % Foundation costs
            if (obj.depth < -60)
                depth = -60;
            else
                depth = obj.depth;
            end
            pile_length = (1+costsParams.embLperDepth) * abs(depth);
            transition_piece = costsParams.transitionLength;
            pile_diameter = sqrt(obj.rating / 10e6) * 9;
            steel_thickness = pile_diameter / costsParams.pileDperThickness * 1000 * (1+costsParams.thickIncr)^((abs(depth) - 30)/5);
            weight = (pile_length + transition_piece) * pile_diameter * pi * steel_thickness/1000 * costsParams.steelWeight;
            steel_costs = weight * costsParams.steelCosts / 1e6;
            hardw_min_steel = costsParams.foundMinSteel * weight / 1000;
            total_no_instal = steel_costs + hardw_min_steel;
            
            % Installation costs
            turbInstall = costsParams.turbInstallation * obj.rating/1e6;
            foundInstall = (costsParams.foundInstallation/1000 * weight);
        
            cost = turbEq + warranty + total_no_instal + turbInstall + foundInstall;
            obj.CAPEX = cost;
        end
    end
end

