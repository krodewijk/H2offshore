classdef Windfarm < handle
    % Group of turbines in predefined grid
    
    properties
        ID;
        dimensions; % Size of farm in intrinsic coordinates
        H2Interconnect; % true --> pipes interconnect, false --> cables interconnect
        centralSubPresent;
        
        % All coordinates of farm
        xIntrin;
        yIntrin;
        
        % x,y coordinates of left bottom
        start_x;
        start_y;
        
        connectionNodes = [];
        
        subPlatform;
        bbConnectionPlatformPresent;
        bbPlatform;
        
        turbines = Turbine.empty; % array of included turbines
        cables = Cable.empty; % turbine cables
        outCables = Cable.empty; % Cable(s) to main hub
        pipes = Pipe.empty; % turbine pipes
        outPipes = Pipe.empty;
        ratedPower = 0; % sum of the rated powers of all connected turbines
        inputPower = 0; % sum of power of all incoming powers via cables/pipes from turbines (year average)
        outputPower = 0; % power output from substation to central hub (year average)
        maxPower; % max allowed connected power
        
        % Turbine cable ratings
        turb2subCableVRating;
        turb2subCableIRating;
        turb2subCableFreq;
        turb2subCableA;
        turb2subCableCap;
        
        % Turbine pipe ratings
        turb2subPipePressure = 30; % bar
        turb2subPipeRadius = 0.25; % m
        turb2subPipeGasFlow; % m3/day
        turb2subPipeEff;
        turb2subInT;
        turb2suboutT;
    end
    
    methods
        function obj = Windfarm(grid, start_x, start_y, central_sub)
            global property;

            obj.start_x = start_x;
            obj.start_y = start_y;
            obj.dimensions = property.farmDim;
            obj.centralSubPresent = central_sub;
            obj.bbConnectionPlatformPresent = property.bbConnectionPlatform;
            
            obj.turb2subCableVRating = property.turb2subCableVRating; % kV 
            obj.turb2subCableIRating = property.turb2subCableIRating; % A
            obj.turb2subCableFreq = property.turb2subCableFreq; % Hz
            obj.turb2subCableA = property.turb2subCableA; % mm2
            obj.turb2subCableCap = property.turb2subCableCap; % uF / km
            
            obj.turb2subPipePressure = property.turb2subPipePressure; % bar 
            obj.turb2subPipeRadius = property.turb2subPipeRadius; % m 
            obj.turb2subPipeGasFlow = property.turb2subPipeGasFlow; % m3/day
            obj.turb2subPipeEff = property.turb2subPipeEff; % 
            obj.turb2subInT = property.turb2subInT; % degC
            obj.turb2suboutT = property.turb2suboutT; % degC
            
            obj.maxPower = property.farmPower * 1e9;
            
            if property.scenario == "H2inTurb"
                obj.H2Interconnect = true;
            else
                obj.H2Interconnect = false;
            end
            
            % Place subplatform (if necessary)
            if obj.centralSubPresent
                obj.placePlatforfm(grid);
            end
            
            % Place turbines
            obj.populate(grid);
                        
            % Connect turbines (with cables or pipes)
            if obj.H2Interconnect
                obj.createPipes(grid);
            else
                if obj.centralSubPresent
                    obj.createDirectCables(grid);
                else
                    obj.createSharedCables(grid);
                end
            end
            
            obj.calculatePower();
        end
        
        function placePlatforfm(obj, grid)
            % Calculate centre
            xCentre = round(obj.start_x + 0.5*obj.dimensions(1));
            yCentre = round(obj.start_y + 0.5*obj.dimensions(2));
            
            % Place platform in middle if possible
            if(grid.mask(xCentre, yCentre))
                % Is possible
                plat = Platform(grid, xCentre, yCentre);
                obj.subPlatform = plat;
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
                                obj.subPlatform = plat;
                                found = true;
                                break;
                            end
                        end
                        if found
                            break;
                        end
                    end
                    radius = radius+1;
                    if (radius > round(0.5*obj.dimensions(1)))
                        error("Couldnt find suitable middle for windfarm");
                    end
                end          
            end
        end
        
        function populate(obj, grid)
            % Create initial coordinates
            xArr = obj.start_x:(obj.start_x + obj.dimensions(1)-1);
            yArr = obj.start_y:(obj.start_y + obj.dimensions(2)-1);
            
            % Make sure coordiantes are within bound
            dim = size(grid.lat);
            maxX = dim(1);
            maxY = dim(2);
            % truncate arrays
            xArr(xArr > maxX) = [];
            yArr(yArr > maxY) = [];
            
            [obj.xIntrin, obj.yIntrin] = meshgrid(xArr, yArr);
            
            turb_ind = 1;
            % Check if coordinates intersect with grid mask and populate
            % with turbines
            for row = obj.start_x:xArr(end)
                % if row is odd increase y
                if (mod(row, 2) == 1)
                    for col = obj.start_y:yArr(end)
                        if(grid.mask(row,col) && (obj.inputPower < obj.maxPower || obj.maxPower <= 0))
                            % Possible turbine location
                            obj.turbines(turb_ind) = Turbine(grid, row, col);
                            % Flip mask pixel
                            grid.mask(row,col) = false;
                            turb_ind = turb_ind + 1;
                        end
                    end
                else
                    % Row is even
                    for col = flip(obj.start_y:yArr(end))
                        if(grid.mask(row,col) && (obj.inputPower < obj.maxPower || obj.maxPower <= 0))
                            obj.turbines(turb_ind) = Turbine(grid, row, col);
                            % Flip mask pixel
                            grid.mask(row,col) = false;
                            turb_ind = turb_ind + 1;
                        end
                    end
                end
            end
        end
        
        function createDirectCables(obj, grid)
            obj.cables = Cable.empty;
            obj.connectionNodes = [];
            
            % Go over the array with turbines
            numTurbs = numel(obj.turbines);
            
            for i = 1:numTurbs
                % Create cable
                cab = Cable(obj.turb2subCableVRating, obj.turb2subCableIRating, obj.turb2subCableFreq, obj.turb2subCableA, obj.turb2subCableCap);

                % Add turbine to cable
                cab.add_turbine(grid, obj.turbines(i));
                cab.update_connections(grid);
                
                % Connect to subplatform
                cab.add_connection(grid, obj.subPlatform.node);
                
                % Save to array
                obj.cables(end+1) = cab;
            end
        end
        
        function createSharedCables(obj, grid)
            obj.cables = Cable.empty;
            obj.connectionNodes = [];
            % Go over the array with turbines
            % It is expected the order of the array is in nice order.
            numTurbs = max(size(obj.turbines));
            cab = Cable(obj.turb2subCableVRating, obj.turb2subCableIRating, obj.turb2subCableFreq, obj.turb2subCableA, obj.turb2subCableCap);
            for i = 1:numTurbs
                % Add turbine to cable
                if ~cab.add_turbine(grid, obj.turbines(i))
                    % No space anymore, save current and create new cable
                    cab.update_connections(grid);
                    obj.cables(end+1) = cab;
                    obj.connectionNodes(1,end+1) = cab.startNode;
                    obj.connectionNodes(2,end) = cab.endNode;
                    cab = Cable(obj.turb2subCableVRating, obj.turb2subCableIRating, obj.turb2subCableFreq, obj.turb2subCableA, obj.turb2subCableCap);
                    cab.add_turbine(grid, obj.turbines(i));
                end
            end
            % Save last cable also...
            cab.update_connections(grid);
            obj.cables(end+1) = cab;
        end
        
        function createPipes(obj, grid)
            % Pipes are shared
            obj.pipes = Pipe.empty;
            obj.connectionNodes = [];
            % Go over the array with turbines
            % It is expected the order of the array is in nice order.
            numTurbs = max(size(obj.turbines));
            pipe = Pipe(obj.turb2subPipeRadius, obj.turb2subPipePressure, obj.turb2subPipeEff, obj.turb2subInT, obj.turb2suboutT);
            for i = 1:numTurbs
                % Add turbine to cable
                if ~pipe.add_turbine(grid, obj.turbines(i))
                    pipe.update_connections(grid);
                    % No space anymore, save current and create new cable
                    obj.pipes(end+1) = pipe;
                    obj.connectionNodes(1,end+1) = pipe.startNode;
                    obj.connectionNodes(2,end) = pipe.endNode;
                    pipe = Pipe(obj.turb2subPipeRadius, obj.turb2subPipePressure, obj.turb2subPipeEff, obj.turb2subInT, obj.turb2suboutT);
                    pipe.add_turbine(grid, obj.turbines(i));
                end
            end
            % Save last cable also...
            pipe.update_connections(grid);
            obj.pipes(end+1) = pipe;
        end
        
        function calculatePower(obj)
            numTurbs = numel(obj.turbines);
            powerSum = 0;
            
            % Sum all the powers of the turbines
            for i = 1:numTurbs
                powerSum = powerSum + obj.turbines(i).rating;
            end
            
            obj.ratedPower = powerSum;
            
            % Sum all the output power of connected cables/pipes
            powerSum = 0;
            if obj.H2Interconnect
                % Sum all pipe powers from turbines
                for i = 1:numel(obj.pipes)
                    powerSum = powerSum + obj.pipes(i).outputPower;
                end
            else
                % Sum all cable powers from turbines
                for i = 1:numel(obj.cables)
                    powerSum = powerSum + obj.cables(i).outputPower;
                end
            end
            
            obj.inputPower = powerSum;
            
            % assume no losses in windfarm
            obj.outputPower = obj.inputPower;
        end
        
        function connect2hub(obj, grid, hubNode)
            if ~obj.centralSubPresent
                % Continue the already laid individual cables
                if obj.H2Interconnect
                    % Connect the pipes
                    numPipes = max(size(obj.pipes));

                    for i = 1:numPipes
                        % Connected turbines
                        turbines = obj.pipes(i).connected.turbines;
                        
                        % Find closest turbine
                        turbCoords = [turbines.xIntrin; turbines.yIntrin]; % 1st row x coord, 2nd row y coord
                        
                        % calculate distanxes to hubNode
                        [hubX, hubY] = grid.node2intrin(hubNode);
                        dists = zeros(1,numel(turbines));
                        
                        for turb = 1:numel(turbines)
                            dists(turb) = norm([hubX; hubY] - turbCoords(:,turb));
                        end
                        
                        [~, id] = min(dists);
                        
                        closestCoords = turbCoords(:,id);
                        closestNode = grid.intrin2node(closestCoords(1), closestCoords(2));
                        
                        new_pipe = Pipe(obj.turb2subPipeRadius, obj.turb2subPipePressure, obj.turb2subPipeEff, obj.turb2subInT, obj.turb2suboutT);
                        
                        new_pipe.add_connection(grid, closestNode);
                        new_pipe.add_connection(grid, hubNode);
                        
                        obj.pipes(end+1) = new_pipe;
                    end
                else
                    % Connect the cables
                    numCables = max(size(obj.cables));

                    for i = 1:numCables
                        obj.cables(i).add_connection(grid, hubNode)
                    end
                end
            else
                % Make cables from the central subplatform to the hub
                
                
            end
        end
        
        function connect2backbone(obj, grid, threshold)
            disp(['Threshold: ', num2str(threshold)]);
            obj.outCables = Cable.empty;
            obj.outPipes = Pipe.empty;
            
            % Find closest backbone pipe
            if obj.centralSubPresent
                hubCoords = [grid.X(obj.subPlatform.xIntrin, obj.subPlatform.yIntrin), grid.Y(obj.subPlatform.xIntrin, obj.subPlatform.yIntrin), 0];
                [hub_lat, hub_lon] = grid.intrin2geo(obj.subPlatform.xIntrin, obj.subPlatform.yIntrin);
            else
                hubXIntrin = obj.start_x + 0.5 * obj.dimensions(1);
                hubYIntrin = obj.start_y + 0.5 * obj.dimensions(2);
                
                hubCoords = [grid.X(hubXIntrin, hubYIntrin), grid.Y(hubXIntrin, hubYIntrin), 0];
                [hub_lat, hub_lon] = grid.intrin2geo(hubXIntrin, hubYIntrin);
            end
                
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
                obj.connect2backbone(grid, threshold + 10);
                return;
            end
            
            end_coords = hubCoords' + shortestVec;
            
            % Convert connection point back to intrinsic values
            [lat, lon] = projinv(grid.projection, end_coords(1), end_coords(2));
            [xIntrin, yIntrin] = grid.geo2intrin(lat, lon);
            connectionNode = grid.intrin2node(xIntrin, yIntrin);
            
            obj.connect2hub(grid, connectionNode);
                
            % Add a platform at the connection point with the backbone
            if obj.bbConnectionPlatformPresent
                obj.bbPlatform = Platform(grid, xIntrin,yIntrin);
            end
        end
        
        function plot(obj, grid)
            %mapaxis = axesm('lambertstd', 'MapLatLimit', grid.latlim, 'MapLonLimit', grid.lonlim);
            % Get limits of map
            xlim = [obj.start_x, obj.start_x + obj.dimensions(1)];
            ylim = [obj.start_y, obj.start_y + obj.dimensions(2)];
            
            numTurbs = max(size(obj.turbines));
            numCables = max(size(obj.cables));
            numPipes = numel(obj.pipes);
            % Plot turbines
            for i = 1:numTurbs
                xCoord = grid.X(obj.turbines(i).xIntrin, obj.turbines(i).yIntrin);
                yCoord = grid.Y(obj.turbines(i).xIntrin, obj.turbines(i).yIntrin);
                
                plot(xCoord, yCoord, 'rx', 'MarkerSize', 20, 'LineWidth',4);
                hold on;
            end
            
            % Plot cables
            for i = 1:numCables
                if any(obj.cables(i).xIntrin)
                    lenCable = max(size(obj.cables(i).xIntrin));
                    xCoord = [];
                    yCoord = [];
                    for k = 1:lenCable
                        xCoord(end+1) = grid.X(obj.cables(i).xIntrin(k), obj.cables(i).yIntrin(k));
                        yCoord(end+1) = grid.Y(obj.cables(i).xIntrin(k), obj.cables(i).yIntrin(k));
                    end
                    plot(xCoord, yCoord,'y', 'LineWidth',3);
                    hold on;
                end
            end
            
            % Plot pipes
            for i = 1:numPipes
                if any(obj.pipes(i).xIntrin)
                    lenPipe = max(size(obj.pipes(i).xIntrin));
                    xCoord = [];
                    yCoord = [];
                    for k = 1:lenPipe
                        xCoord(end+1) = grid.X(obj.pipes(i).xIntrin(k), obj.pipes(i).yIntrin(k));
                        yCoord(end+1) = grid.Y(obj.pipes(i).xIntrin(k), obj.pipes(i).yIntrin(k));
                    end
                    plot(xCoord, yCoord,'c', 'LineWidth',3);
                    hold on;
                end
            end
            
            % Plot central subplatform
            if obj.centralSubPresent
                plot(grid.X(obj.subPlatform.xIntrin, obj.subPlatform.yIntrin), grid.Y(obj.subPlatform.xIntrin, obj.subPlatform.yIntrin), "sb", 'MarkerSize', 10, 'LineWidth', 7);
            end
            
            if numel(obj.bbPlatform) > 0
                plot(grid.X(obj.bbPlatform.xIntrin, obj.bbPlatform.yIntrin), grid.Y(obj.bbPlatform.xIntrin, obj.bbPlatform.yIntrin), "sb", 'MarkerSize', 10, 'LineWidth', 7);
            end
            
            % Plot bounding box
            xBox = [obj.start_x, obj.start_x, obj.start_x + obj.dimensions(1)-1, obj.start_x + obj.dimensions(1)-1, obj.start_x];
            yBox = [obj.start_y, obj.start_y + obj.dimensions(2)-1, obj.start_y + obj.dimensions(2)-1, obj.start_y, obj.start_y];
            
            for i = 1:numel(xBox)
                xBoxMap(i) = grid.X(xBox(i), yBox(i));
                yBoxMap(i) = grid.Y(xBox(i), yBox(i));
            end
            
            plot(xBoxMap, yBoxMap, '--', 'LineWidth',2);
            hold off;
        end
        
        function plotIntrin(obj)
            % Get limits of map
            xlim = [obj.start_x, obj.start_x + obj.dimensions(1)];
            ylim = [obj.start_y, obj.start_y + obj.dimensions(2)];
            
            numTurbs = max(size(obj.turbines));
            numCables = max(size(obj.cables));
            % Plot turbines
            for i = 1:numTurbs
                plot(obj.turbines(i).xIntrin, obj.turbines(i).yIntrin, 'rx', 'MarkerSize', 20, 'LineWidth',4);
                hold on;
            end
            
            % Plot cables
            for i = 1:numCables
                plot(obj.cables(i).xIntrin, obj.cables(i).yIntrin,'c', 'LineWidth',5);
                hold on;
            end
            hold off;
        end
    end
end

