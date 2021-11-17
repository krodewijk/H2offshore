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
        wind_hubheight; % m/s
        windVar;
        windDist;
        windDistBins;
        power; % W (actually generated power)
        hubHeight; % m
        H2out; % true for inturbine electrolyser, false for electric output
        H2Eff;
        floating;

        xIntrin; % with reference to defined grid 
        yIntrin;
        
        costOffshore = nan;
        costOnLoc = nan;
        
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
            
            % Save mean wind speed at hub height for future reference
            obj.wind_hubheight = hist * x;
            
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
        
        function cost = calculateCostFromExcel(obj, costTables, distance)
            % Calculates the costs corresponding to one wind turbine from
            % the data supplied in the excel sheet. The transportation to
            % shore is also taken into account.
            
            % determine if the turbine is floating or fixed
            breakevenPoint = interp1(costTables.distances.Variables, costTables.breakevens.Variables, distance);
            if abs(obj.depth) > breakevenPoint
                % floating
                obj.floating = true;
            else
                % fixed
                obj.floating = false;
            end

            % find the locations in the cost table and interpolate
            if obj.H2out
                % hydrogen output
                if obj.floating
                    % and floating
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.float.h2onshore.Variables, double(distance), double(abs(obj.depth)));
                else
                    % and fixed
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.fixed.h2onshore.Variables, double(distance), double(abs(obj.depth)));
                end
            else
                % electricity output
                if obj.floating
                    % and floating
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.float.grid.Variables, double(distance), double(abs(obj.depth)));
                else
                    % and fixed
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.fixed.grid.Variables, double(distance), double(abs(obj.depth)));
                end
            end
            obj.costOffshore = cost;
        end

        function cost = calculateCostFromExcelOnLoc(obj, costTables, distance)
            % Calculates the costs corresponding to one wind turbine from
            % the data supplied in the excel sheet. The transportation to
            % shore is not taken into account
            
            % determine if the turbine is floating or fixed
            breakevenPoint = interp1(costTables.distances.Variables, costTables.breakevens.Variables, distance);
            if abs(obj.depth) > breakevenPoint
                % floating
                obj.floating = true;
            else
                % fixed
                obj.floating = false;
            end

            % find the locations in the cost table and interpolate
            if obj.H2out
                % hydrogen output
                if obj.floating
                    % and floating
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.float.h2offshore.Variables, double(distance), double(abs(obj.depth)));
                else
                    % and fixed
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.fixed.h2offshore.Variables, double(distance), double(abs(obj.depth)));
                end
            else
                % electricity output
                if obj.floating
                    % and floating
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.float.turbine.Variables, double(distance), double(abs(obj.depth)));
                else
                    % and fixed
                    cost = interp2(costTables.distances.Variables, costTables.depths.Variables, costTables.fixed.turbine.Variables, double(distance), double(abs(obj.depth)));
                end
            end
            obj.costOnLoc = cost;
        end
    end
end

