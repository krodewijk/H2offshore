classdef Cluster < handle
    % Windfarms and hub
    
    properties
         scenario;
         mappingOption;
         minTurbinesPerFarm;
         farmDim;
         dimension;
         maxPower; % only used if mappingOption = "power"
         
         inputPower; % Power coming in from cables
         outputPower; % Power going out from central hub to shore
         
         % Connected objects
         platform = Platform.empty;
         windfarms = Windfarm.empty;
         cables2shore = Cable.empty;
         pipes2shore = Pipe.empty;
         
         % Cable properties
         hub2shoreCableVRating; % kV (voltage rating)
         hub2shoreCableIRating; % A (current rating)
         hub2shoreCableFreq; % if cable is AC or DC (true = AC, false = DC)
         hub2shoreCableCross; % cross section
         hub2shoreCableCap; % capacitance / km
         
         % Pipe properties
         sub2mainPipePressure;
         sub2mainPipeRadius;
         sub2mainPipeEff;
         sub2mainInT;
         sub2mainoutT;
         
         xlim;
         ylim;
    end
    
    methods
        function obj = Cluster()
            % Constructor
            global property;
            obj.farmDim = property.farmDim;
            obj.minTurbinesPerFarm = property.farmMinTurbs;
            
            obj.hub2shoreCableVRating = property.hub2shoreCableVRating;
            obj.hub2shoreCableIRating = property.hub2shoreCableIRating;
            obj.hub2shoreCableFreq = property.sub2hubCableFreq;
            obj.hub2shoreCableCross = property.sub2hubCableA;
            obj.hub2shoreCableCap = property.hub2shoreCableCap;
            
            obj.scenario = property.scenario;
            
            obj.sub2mainPipePressure = property.sub2mainPipePressure; % bar (pressure rating)
            obj.sub2mainPipeRadius = property.sub2mainPipeRadius; % m (pipe radius)
            obj.sub2mainPipeEff = property.sub2mainPipeEff;
            obj.sub2mainInT = property.sub2mainInT;
            obj.sub2mainoutT = property.sub2mainoutT;
            
            obj.mappingOption = property.clusterMappingOption;
            obj.dimension = property.clusterDim;
            obj.maxPower = property.clusterPower * 1e9;
        end
        
        function populate(obj, grid, x, y)
            grid.resetMask();
            if obj.mappingOption == "size"
                % x,y now used as cluster starting coordinates
                % (Bottom-left)
                xlim = [x, x + obj.dimension(1)];
                ylim = [y, y + obj.dimension(2)];
                obj.populateRegion(grid, xlim, ylim);
            elseif obj.mappingOption == "power"
                % x,y now used as platform coordinates
                obj.populateAround(grid, x, y);
            else
                error("Given mapping option not valid. (check parameters)");
            end
        end
        
        function populateAround(obj, grid, xPlatform, yPlatform)
            % Create cluster from scratch around platform until max power
            % is reached
            obj.platform = Platform.empty;
            obj.windfarms = Windfarm.empty;
            obj.cables2shore = Cable.empty;
            obj.pipes2shore = Pipe.empty;
            
            % Calculate limits of region
            obj.xlim = [xPlatform - round(obj.dimension(1)/2), xPlatform + round(obj.dimension(1)/2)];
            obj.ylim = [yPlatform - round(obj.dimension(2)/2), yPlatform + round(obj.dimension(2)/2)];
        
            % Check if platform coords is valid
            if (grid.mask(xPlatform, yPlatform))
                % Coords ok, place platform
                obj.platform = Platform(grid, xPlatform, yPlatform);
            else
                error("Given platform coordinates are not allowed!");
            end
            
            radius = 1;
            powerLimitReached = false;
            
            % Place first windfarm at
            newFarm = Windfarm(grid, xPlatform, yPlatform, false);
            % Check if farm has more than property.minTurbs
            % turbines
            if max(size(newFarm.turbines)) >= obj.minTurbinesPerFarm
                obj.windfarms(end+1) = newFarm;
                obj.windfarms(end).connect2hub(grid, obj.platform.node);
            end
            
            % Check if power limit is reached
            obj.calculatePower();
            if obj.inputPower > obj.maxPower
                powerLimitReached = true;
            end
            
            % Place windfarms around platform
            while (radius < floor(min(obj.dimension)/2)) && ~powerLimitReached
                % Start in top right corner
                curOffset = [radius, radius];
                % For each edge (of 4 edges)
                for edge = 1:4
                    if edge == 1
                        % right edge (going down)
                        step = [0, -1];
                    elseif edge == 2
                        % bottom edge (going left)
                        step = [-1, 0];
                    elseif edge == 3
                        % left edge (going up)
                        step = [0, 1];
                    elseif edge == 4
                        % top edge (going right)
                        step = [1, 0];
                    end
                    for tim = 1:(2*radius)
                        curOffset = curOffset + step;
                        curCoords = curOffset .* obj.farmDim + [xPlatform, yPlatform];
                        % Place farm
                        newFarm = Windfarm(grid, curCoords(1), curCoords(2), false);
                        % Check if farm has more than property.minTurbs
                        % turbines
                        if max(size(newFarm.turbines)) >= obj.minTurbinesPerFarm
                            obj.windfarms(end+1) = newFarm;
                            obj.windfarms(end).connect2hub(grid, obj.platform.node);
                        end
                        % Check if power limit is reached
                        obj.calculatePower();
                        if obj.inputPower > obj.maxPower
                            powerLimitReached = true;
                            % Remove last windfarm
                            %obj.windfarms(end) = [];
                            break;
                        end
                    end
                    if powerLimitReached
                        break;
                    end
                end
                radius = radius + 1;
            end
            obj.calculatePower();
        end
        
        function populateRegion(obj, grid, xlim, ylim)
            % Create cluster from scratch in assigned region
            obj.platform = Platform.empty;
            obj.windfarms = Windfarm.empty;
            obj.cables2shore = Cable.empty;
            obj.pipes2shore = Pipe.empty;
            obj.xlim = xlim;
            obj.ylim = ylim;
            
            % Calculate middle of cluster
            xCentre = min(xlim) + round((max(xlim) - min(xlim)) / 2);
            yCentre = min(ylim) + round((max(ylim) - min(ylim)) / 2);
            
            % Place platform in middle if possible
            if(grid.mask(xCentre, yCentre))
                % Is possible
                plat = Platform(grid, xCentre, yCentre);
                obj.platform = plat;
            else
                % Need to find other place...
                found = false;
                radius = 1;
                while ~found 
                    curOffset = [radius, radius]; % start at top right
                    % For each edge (of 4 edges)
                    for edge = 1:4
                        if edge == 1
                            % right edge (going down)
                            step = [0, -1];
                        elseif edge == 2
                            % bottom edge (going left)
                            step = [-1, 0];
                        elseif edge == 3
                            % left edge (going up)
                            step = [0, 1];
                        elseif edge == 4
                            % top edge (going right)
                            step = [1, 0];
                        end
                        for tim = 1:(2*radius)
                            curOffset = curOffset + step;
                            curCoords = curOffset + [xCentre, yCentre];
                            % Place platform in middle if possible
                            if(grid.mask(curCoords(1), curCoords(2)))
                                % Is possible
                                plat = Platform(grid, curCoords(1), curCoords(2));
                                obj.platform = plat;
                                found = true;
                                break;
                            end
                        end
                        if found
                            break;
                        end
                    end
                    radius = radius+1;
                    if (radius > round(max(xlim)/2))
                        error("Couldnt find suitable middle for cluster");
                    end
                end          
            end
            
            fprintf(1,'\nPlatform placed at [%i, %i]',obj.platform.xIntrin, obj.platform.yIntrin);          
            
            % Platform placed!
            % Time to place windfarms
            xFarms = floor(abs(max(xlim) - min(xlim)) / obj.farmDim(1));
            yFarms = floor(abs(max(ylim) - min(ylim)) / obj.farmDim(2));
            
            xStep = obj.farmDim(1);
            yStep = obj.farmDim(2);
            
            startX = min(xlim):xStep:(max(xlim)-obj.farmDim(1));
            startY = min(ylim):yStep:(max(ylim)-obj.farmDim(2));
            
            fprintf(1,'\nCreating windfarms\nProgress:    ');
            for x = 1:xFarms
                for y = 1:yFarms
                    % Create farm
                    newFarm = Windfarm(grid, startX(x), startY(y), false);
                    % Check if farm has more than property.minTurbs
                    % turbines
                    if max(size(newFarm.turbines)) >= obj.minTurbinesPerFarm
                        obj.windfarms(end+1) = newFarm;
                        obj.windfarms(end).connect2hub(grid, obj.platform.node);
                    end
                end
                fprintf(1,'\b\b\b\b%3.0f%%',(x/xFarms*100));
            end
            fprintf(1,'\n');
            obj.calculatePower();
        end
        
        function calculatePower(obj)
            numFarms = numel(obj.windfarms);
            powerSum = 0;
            
            % Sum all the powers coming in from the cables
            for i = 1:numFarms
                powerSum = powerSum + obj.windfarms(i).outputPower;
            end
            
            obj.inputPower = powerSum;
            
            % Calculate power going out
            obj.outputPower = obj.platform.calculate_power(obj.inputPower);
        end
        
        function connect2shore(obj, grid, option)
            obj.cables2shore = Cable.empty;
            obj.pipes2shore = Pipe.empty;
            possibleOptions = ["closest", "NL", "DK", "DE", "NOR", "UK", "BE", "FR"];
            if ~any(possibleOptions == option)
                % given option not possible
                error("Given option not possible.");
            end
            
            % Go to closest connection point
            numPoints = numel(grid.shapes.connectionPoints);
            dists = nan(numPoints, 1);
            for point = 1:numPoints
                % Skip current point is node == NaN
                if isnan(grid.shapes.connectionPoints(point).node)
                    continue;
                end
                if option == "closest"
                    % Calculate distance to every onshore point
                    dists(point) = distances(grid.graph, obj.platform.node, grid.shapes.connectionPoints(point).node);
                elseif grid.shapes.connectionPoints(point).Country == option
                    dists(point) = distances(grid.graph, obj.platform.node, grid.shapes.connectionPoints(point).node);
                else
                    dists(point) = NaN;
                end
            end
            [~, ind] = min(dists); % calc which point is closest
            closestPoint = grid.shapes.connectionPoints(ind);
            
            if obj.scenario == "fullElectric"
                % Create cable
                cable = Cable(obj.hub2shoreCableVRating, obj.hub2shoreCableIRating, obj.hub2shoreCableFreq, obj.hub2shoreCableCross, obj.hub2shoreCableCap);
                % Add connections
                cable.add_connection(grid, obj.platform.node);
                cable.add_connection(grid, closestPoint.node);

                numCables = ceil(obj.outputPower / cable.power_rating);
                for i = 1:numCables
                    obj.cables2shore(end+1) = cable;
                end
            else
                % Create pipe
                pipe = Pipe(obj.sub2mainPipeRadius, obj.hub2shorePipePressure, obj.hub2shorePipeGasV);
                % Add connections
                pipe.add_connection(grid, obj.platform.node);
                pipe.add_connection(grid, closestPoint.node);

                numPipes = ceil(obj.outputPower / pipe.power_rating);
                for i = 1:numPipes
                    obj.pipes2shore(end+1) = pipe;
                end
            end
        end
        
        function connect2backbone(obj, grid, threshold)
            disp(['Threshold: ', threshold]);
            obj.cables2shore = Cable.empty;
            obj.pipes2shore = Pipe.empty;
            
            % Find closest backbone pipe
            hubCoords = [grid.X(obj.platform.xIntrin, obj.platform.yIntrin), grid.Y(obj.platform.xIntrin, obj.platform.yIntrin), 0];
            [hub_lat, hub_lon] = grid.intrin2geo(obj.platform.xIntrin, obj.platform.yIntrin);
            
            numPipes = numel(grid.shapes.pipelines);
            
            shortestDist = nan;
            shortestVec = nan;
            
            % Find the shortest vector from the hub to the pipeline
            % sections
            fprintf(1,'\n Connecting to backbone:    ');
            for pipe = 1:numPipes
                % For each pipe in the backbone
                pipeline = grid.shapes.pipelines(pipe);
                numSections = numel(pipeline.Lat);
                
                for section = 1:numSections-1
                    % For each section of the pipe
                    if isnan(pipeline.Lat(section)) || isnan(pipeline.Lat(section+1))
                        continue;
                    end
                    
                    % Store the current and next coordinates
                    cur_lat = pipeline.Lat(section);
                    next_lat = pipeline.Lat(section+1);
                    cur_lon = pipeline.Lon(section);
                    next_lon = pipeline.Lon(section+1);
                                                            
                    % Check if section is sufficiently close for the
                    % calculation (takes a longg time otherwise...)
                    sphere_dist = (distance(hub_lat, hub_lon, cur_lat, cur_lon) / 360)*2*pi*6371;
                    if  sphere_dist > threshold
                        % Point on pipeline is too far away
                        continue;
                    end
                    
                    % Calculate the vector coords in km
                    [x,y] = projfwd(grid.projection, cur_lat, cur_lon);
                    cur = [x, y, 0];
                    [x,y] = projfwd(grid.projection, next_lat, next_lon);
                    next = [x, y, 0];
                    
                    
                    % Calculate triangle vectors AB and AC (A,B = pipe
                    % points. C = hub vector)
                    AB = (next - cur)';
                    BA = -AB;
                    AC = (hubCoords - cur)';
                    BC = (hubCoords - next)';
                    
                    % Calculate the orthogonal distance d between the hub and
                    % the pipe section
                    if dot(AB, AC) < 0
                        % Point lies before A
                        % Shortest vector is -AC
                        % vec is from hub(C) to the pipeline
                        vec = -AC;
                        d = norm(AC);
                    elseif dot(BA, BC) < 0
                        % Point lies after B
                        % Shortest vector is BC
                        vec = -BC;
                        d = norm(vec);
                    else
                        % Point lies between A and B
                        d = norm(cross(AB, AC)) / norm(AB);
                        AE_mag = sqrt(norm(AC)^2 - d^2);
                        AE = AB ./ norm(AB) .* AE_mag;
                        vec = AE - AC;
                    end
                    
                    if d < shortestDist || isnan(shortestDist)
                        % Shortest vector till now, store it
                        shortestDist = d;
                        shortestVec = vec;
                    end
                end
                fprintf(1,'\b\b\b\b%3.0f%%',(pipe/numPipes*100));
            end
            
            % Check if shortest vector is found otherwise try again with a
            % higher threshold
            if isnan(shortestDist)
                obj.connect2backbone(grid, threshold + 50);
                return;
            end
            
            end_coords = hubCoords + shortestVec;
            
            % Convert connection point back to intrinsic values
            [lat, lon] = projinv(grid.projection, end_coords(1), end_coords(2));
            [xIntrin, yIntrin] = grid.geo2intrin(lat, lon);
            connectionNode = grid.intrin2node(xIntrin, yIntrin);
            
            if obj.scenario == "fullElectric"
                % Create cable
                cable = Cable(obj.hub2shoreCableVRating, obj.hub2shoreCableIRating, obj.hub2shoreCableAC);
                % Add connections
                cable.add_connection(grid, obj.platform.node);
                cable.add_connection(grid, connectionNode);

                numCables = ceil(obj.outputPower / cable.power_rating);
                for i = 1:numCables
                    obj.cables2shore(end+1) = cable;
                end
            else
                % Create pipe
                pipe = Pipe(obj.sub2mainPipeRadius, obj.sub2mainPipePressure, obj.sub2mainPipeGasV);
                % Add connections
                pipe.add_connection(grid, obj.platform.node);
                pipe.add_connection(grid, connectionNode);

                numPipes = ceil(obj.outputPower / pipe.power_rating);
                for i = 1:numPipes
                    obj.pipes2shore(end+1) = pipe;
                end
            end
        end
        
        function plot(obj, grid)
            coastLines = load("coastlines");
            cities = shaperead('worldcities', 'UseGeoCoords', true);
            citiesLat = zeros(max(size(cities),1));
            citiesLon = zeros(max(size(cities)),1);
            
            % Convert shape format to useable array
            for i = 1:max(size(cities))
                citiesLat(i) = cities(i).Lat;
                citiesLon(i) = cities(i).Lon;
            end
            
            % Plot windfarms
            numFarms = max(size(obj.windfarms));
            for i = 1:numFarms
                obj.windfarms(i).plot(grid);
                hold on;
            end
            
            % Plot cable2shore
            numCables = numel(obj.cables2shore);
            for i = 1:numCables
                lenCable = max(size(obj.cables2shore(i).xIntrin));
                xCoord = [];
                yCoord = [];
                for k = 1:lenCable
                    xCoord(end+1) = grid.X(obj.cables2shore(i).xIntrin(k), obj.cables2shore(i).yIntrin(k));
                    yCoord(end+1) = grid.Y(obj.cables2shore(i).xIntrin(k), obj.cables2shore(i).yIntrin(k));
                end
                plot(xCoord, yCoord,'y', 'LineWidth',5);
                hold on;
            end
            
            % Plot pipes2shore
            numPipes = numel(obj.pipes2shore);
            for i = 1:numPipes
                lenCable = max(size(obj.pipes2shore(i).xIntrin));
                xCoord = [];
                yCoord = [];
                for k = 1:lenCable
                    xCoord(end+1) = grid.X(obj.pipes2shore(i).xIntrin(k), obj.pipes2shore(i).yIntrin(k));
                    yCoord(end+1) = grid.Y(obj.pipes2shore(i).xIntrin(k), obj.pipes2shore(i).yIntrin(k));
                end
                plot(xCoord, yCoord,'b', 'LineWidth',5);
                hold on;
            end
            
            % Plot platform
            plot(grid.X(obj.platform.xIntrin, obj.platform.yIntrin), grid.Y(obj.platform.xIntrin, obj.platform.yIntrin), "sb", 'MarkerSize', 10, 'LineWidth', 7);
            
            % Plot surroundings (coastlines + cities)
            [xProjCoast, yProjCoast] = projfwd(grid.projection, coastLines.coastlat, coastLines.coastlon);
            [xProjCity, yProjCity] = projfwd(grid.projection, citiesLat, citiesLon);
            plot(xProjCoast, yProjCoast);
            plot(xProjCity, yProjCity, 'o', 'MarkerSize', 10, 'LineWidth', 5);
            
            
            
            xlim([-2e5, 7.5e5]);
            ylim([5.4e6, 7.4e6]);
        end
    end
end

