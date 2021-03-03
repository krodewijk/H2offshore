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
                    obj.cables(end+1) = cab;
                    obj.connectionNodes(1,end+1) = cab.startNode;
                    obj.connectionNodes(2,end) = cab.endNode;
                    cab = Cable(obj.turb2subCableVRating, obj.turb2subCableIRating, obj.turb2subCableFreq, obj.turb2subCableA, obj.turb2subCableCap);
                    cab.add_turbine(grid, obj.turbines(i));
                end
            end
            % Save last cable also...
            obj.cables(end+1) = cab;
        end
        
        function createPipes(obj, grid)
            obj.pipes = Pipe.empty;
            obj.connectionNodes = [];
            % Go over the array with turbines
            % It is expected the order of the array is in nice order.
            numTurbs = max(size(obj.turbines));
            pipe = Pipe(obj.turb2hubPipeRadius, obj.turb2hubPipePressure, obj.turb2hubPipeGasV);
            for i = 1:numTurbs
                % Add turbine to cable
                if ~pipe.add_turbine(grid, obj.turbines(i))
                    % No space anymore, save current and create new cable
                    obj.pipes(end+1) = pipe;
                    obj.connectionNodes(1,end+1) = pipe.startNode;
                    obj.connectionNodes(2,end) = pipe.endNode;
                    pipe = Pipe(obj.turb2hubPipeRadius, obj.turb2hubPipePressure, obj.turb2hubPipeGasV);
                    pipe.add_turbine(grid, obj.turbines(i));
                end
            end
            % Save last cable also...
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
                        obj.pipes(i).add_connection(grid, hubNode)
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

