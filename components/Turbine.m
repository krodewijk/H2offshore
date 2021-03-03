classdef Turbine < handle
    % Class decribing a single turbine
    
    properties
        rating; % W (power rating)
        density = 1.225; % kg/m3
        A; % m^2
        max_wind; % m/s
        
        depth; % m
        wind; % m/s
        windVar;
        windDist;
        windDistBins;
        power; % W (actually generated power)
        hubHeight; % m
        H2out; % true for inturbine electrolyser, false for electric output
        
        xIntrin; % with reference to defined grid 
        yIntrin;
        
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
                curve(i) = obj.calculatePower(x(i));
            end
            
            % Multiply random distribution with curve and integrate
            expected = trapz(x, dist.*curve);
            
            obj.power = expected;
        end
        
        function turb_power = calculatePower(obj, wind)
            if (wind < obj.max_wind)
                % Calculate wind power
                wind_power = 0.5 * obj.density * obj.A * wind^3;

                % Assume efficiency of 45%
                turb_power = 0.45 * wind_power;

                % Adjust for rating
                turb_power(turb_power > obj.rating) = obj.rating;

            else
                turb_power = 0;
            end
        end
    end
end

