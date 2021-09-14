classdef Windfarm < handle
    % Group of turbines in predefined grid
    
    properties
        ID;
        dimensions; % Size of farm in intrinsic coordinates
        H2Interconnect; % true --> pipes interconnect, false --> cables interconnect
        centralSubPresent;
        scenario; %'H2inTurb', 'fullElectric' or 'electricToH2'
        onshoreConnection;
        
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
        bbLength = 0;
        bbLines;
        bbNodes;
        bbTransport;
        
        turbines = Turbine.empty; % array of included turbines
        cables = Cable.empty; % turbine cables
        outCables = Cable.empty; % Cable(s) to main hub
        pipes = Pipe.empty; % turbine pipes
        outPipes = Pipe.empty;
        
        capacityFactor = 0;
        
        ratedPower = 0; % sum of the rated powers of all connected turbines
        inputPower = 0; % sum of power of all incoming powers via cables/pipes from turbines (year average)
        outputPower = 0; % power output from substation to central hub (year average)
        shorePower = 0;
        
        inputEnergy = 0; % Expected generated energy over 1 year in GWh
        outputEnergy = 0; % Expected generated energy that reaches backbone over 1 year in GWh
        shoreEnergy = 0; % Expected generated energy that reaches shore over 1 year in GWh
        
        maxPower; % max allowed connected power
        mean_wind; % average wind speed in windpark (m/s)
        mean_depth; % average water depth (m)
        
        transport2shoreLosses;
        CAPEX;
        CAPEXlist;
        OPEX;
        OPEXlist;
        LCOE; % EUR / MWh
        LCOElist;
        LCOH; % EUR / kg
        H2Production; % kg
        
        % Turbine cable ratings
        turb2subCableVRating;
        turb2subCableIRating;
        turb2subCableFreq;
        turb2subCableA;
        turb2subCableCap;
        
        sub2mainCableVRating;
        sub2mainCableIRating;
        sub2mainCableFreq;
        sub2mainCableA;
        sub2mainCableCap;
        
        
        hub2shoreCableVRating;
        hub2shoreCableIRating;
        hub2shoreCableFreq;
        hub2shoreCableA;
        hub2shoreCableCap;
        
        % Turbine pipe ratings
        turb2subPipePressure = 30; % bar
        turb2subPipeRadius = 0.25; % m
        turb2subPipeGasFlow; % m3/day
        turb2subPipeEff;
        turb2subInT;
        turb2suboutT;
        
        sub2mainPipePressure; % bar
        sub2mainPipeRadius; % m
        sub2mainPipeGasFlow; % m3/day
        sub2mainPipeEff;
        sub2mainInT;
        sub2mainoutT;
        
        hub2shorePipePressure; % bar
        hub2shorePipeRadius; % m
        hub2shorePipeEff;
        hub2shoreInT;
        hub2shoreoutT;
    end
    
    methods
        function obj = Windfarm(grid, start_x, start_y, central_sub)
            global property;

            obj.start_x = start_x;
            obj.start_y = start_y;
            obj.dimensions = property.farmDim;
            obj.centralSubPresent = central_sub;
            obj.bbConnectionPlatformPresent = property.bbConnectionPlatform;
            obj.onshoreConnection = property.shoreConnectionPoint;
            
            obj.turb2subCableVRating = property.turb2subCableVRating; % kV 
            obj.turb2subCableIRating = property.turb2subCableIRating; % A
            obj.turb2subCableFreq = property.turb2subCableFreq; % Hz
            obj.turb2subCableA = property.turb2subCableA; % mm2
            obj.turb2subCableCap = property.turb2subCableCap; % uF / km
            
            obj.sub2mainCableVRating = property.sub2hubCableVRating; % kV 
            obj.sub2mainCableIRating = property.sub2hubCableIRating; % A
            obj.sub2mainCableFreq = property.sub2hubCableFreq; % Hz
            obj.sub2mainCableA = property.sub2hubCableA; % mm2
            obj.sub2mainCableCap = property.sub2hubCableCap; % uF / km
            
            obj.hub2shoreCableVRating = property.hub2shoreCableVRating;
            obj.hub2shoreCableIRating = property.hub2shoreCableIRating;
            obj.hub2shoreCableFreq = property.hub2shoreCableFreq;
            obj.hub2shoreCableA = property.hub2shoreCableA;
            obj.hub2shoreCableCap = property.hub2shoreCableCap;
            
            obj.turb2subPipePressure = property.turb2subPipePressure; % bar 
            obj.turb2subPipeRadius = property.turb2subPipeRadius; % m 
            obj.turb2subPipeGasFlow = property.turb2subPipeGasFlow; % m3/day
            obj.turb2subPipeEff = property.turb2subPipeEff; % 
            obj.turb2subInT = property.turb2subInT; % degC
            obj.turb2suboutT = property.turb2suboutT; % degC
            
            obj.sub2mainPipePressure = property.sub2mainPipePressure; % bar 
            obj.sub2mainPipeRadius = property.sub2mainPipeRadius; % m 
            obj.sub2mainPipeGasFlow = property.sub2mainPipeGasFlow; % m3/day
            obj.sub2mainPipeEff = property.sub2mainPipeEff; % 
            obj.sub2mainInT = property.sub2mainInT; % degC
            obj.sub2mainoutT = property.sub2mainoutT; % degC
            
            obj.hub2shorePipePressure = property.main2bbPipePressure; % bar
            obj.hub2shorePipeRadius = property.main2bbPipeRadius; % m
            obj.hub2shorePipeEff = property.main2bbPipeEff;

            obj.scenario = property.scenario;
            
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
            
            % calculate mean wind speed
            obj.calc_mean_wind();
            obj.calc_mean_depth();
        end
        
        function placePlatforfm(obj, grid)
            % Calculate centre
            xCentre = round(obj.start_x + 0.5*obj.dimensions(1));
            yCentre = round(obj.start_y + 0.5*obj.dimensions(2));
            
            % Place platform in middle if possible
            if(grid.mask(xCentre, yCentre))
                % Is possible
                if obj.H2Interconnect
                    plat = Platform(grid, xCentre, yCentre, "H2toH2", "sub");
                else
                    plat = Platform(grid, xCentre, yCentre, "EtoH2", "sub");
                end
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
                                if obj.H2Interconnect
                                    plat = Platform(grid, curCoords(1), curCoords(2), "H2toH2", "sub");
                                else
                                    plat = Platform(grid, curCoords(1), curCoords(2), "EtoE", "sub");
                                end
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
            usedPower = 0;
            % Check if coordinates intersect with grid mask and populate
            % with turbines
            for row = obj.start_x:xArr(end)
                % if row is odd increase y
                if (mod(row, 2) == 1)
                    for col = obj.start_y:yArr(end)
                        if(grid.mask(row,col) && (usedPower < obj.maxPower || obj.maxPower <= 0))
                            % Possible turbine location
                            obj.turbines(turb_ind) = Turbine(grid, row, col);
                            % Add used power
                            usedPower = usedPower + obj.turbines(turb_ind).rating;
                            % Flip mask pixel
                            grid.mask(row,col) = false;
                            turb_ind = turb_ind + 1;
                        end
                    end
                else
                    % Row is even
                    for col = flip(obj.start_y:yArr(end))
                        if(grid.mask(row,col) && (usedPower < obj.maxPower || obj.maxPower <= 0))
                            obj.turbines(turb_ind) = Turbine(grid, row, col);
                            % Add used power
                            usedPower = usedPower + obj.turbines(turb_ind).rating;
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
                        
                        new_pipe.connected.turbines = turbines;
                        new_pipe.power_capacity = obj.pipes(i).power_capacity;
                        new_pipe.add_connection(grid, closestNode);
                        new_pipe.add_connection(grid, hubNode);
                        
                        obj.outPipes(end+1) = new_pipe;
                    end
                else
                    % Connect the cables
                    numCables = max(size(obj.cables));

                    for i = 1:numCables
                        % Connected turbines
                        turbines = obj.cables(i).connected.turbines;
                        
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
                        
                        new_cable = Cable(obj.turb2subCableVRating, obj.turb2subCableIRating, obj.turb2subCableFreq, obj.turb2subCableA, obj.turb2subCableCap);
                            
                        new_cable.connected.turbines = turbines;
                        new_cable.power_capacity = obj.cables(i).power_capacity;
                        new_cable.add_connection(grid, closestNode);
                        new_cable.add_connection(grid, hubNode);
                        
                        obj.outCables(end+1) = new_cable;
                    end
                end
            else
                % Make cables/pipes from the central subplatform to the hub
                if obj.subPlatform.kind == "EtoH2" || obj.subPlatform.kind == "H2toH2"
                    new_pipe = Pipe(obj.sub2mainPipeRadius, obj.sub2mainPipePressure, obj.sub2mainPipeEff, obj.sub2mainInT, obj.sub2mainoutT);
                    % Calculate num of pipes needed to transport all energy
                    numPipes = ceil(obj.inputPower / new_pipe.power_rating);
                    for i = 1:numPipes
                        new_pipe = Pipe(obj.sub2mainPipeRadius, obj.sub2mainPipePressure, obj.sub2mainPipeEff, obj.sub2mainInT, obj.sub2mainoutT);
                        new_pipe.add_connection(grid, obj.subPlatform.node);
                        new_pipe.add_connection(grid, hubNode);

                        obj.outPipes(end+1) = new_pipe;
                    end
                else
                    new_cable = Cable(obj.sub2mainCableVRating, obj.sub2mainCableIRating, obj.sub2mainCableFreq, obj.sub2mainCableA, obj.sub2mainCableCap);
                    % Calculate num of pipes needed to transport all energy
                    numCables = ceil(obj.inputPower / new_cable.power_rating);
                    for i = 1:numCables
                        new_cable = Cable(obj.sub2mainCableVRating, obj.sub2mainCableIRating, obj.sub2mainCableFreq, obj.sub2mainCableA, obj.sub2mainCableCap);
                        new_cable.add_connection(grid, obj.subPlatform.node);
                        new_cable.add_connection(grid, hubNode);

                        obj.outCables(end+1) = new_cable;
                    end
                end
            end
        end
        
        function connect2backbone(obj, grid, threshold)
            %threshold = 50; % km
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
            %fprintf(1,'\n Connecting to backbone:    ');
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
                        chosenPipe = pipe;
                        chosenSection = section;
                    end
                end
                %fprintf(1,'\b\b\b\b%3.0f%%',(pipe/numPipes*100));
            end
            
            % Check if shortest vector is found otherwise try again with a
            % higher threshold
            if isnan(shortestDist)
                obj.connect2backbone(grid, threshold+50);
                return;
            end
            
            end_coords = hubCoords' + shortestVec;
            
            % Convert connection point back to intrinsic values
            [lat, lon] = projinv(grid.projection, end_coords(1), end_coords(2));
            [xIntrin, yIntrin] = grid.geo2intrin(lat, lon);
            connectionNode = grid.intrin2node(xIntrin, yIntrin);
            
            
            if ~isnan(connectionNode)
                obj.connect2hub(grid, connectionNode);
            else
                warning("Can't connect to backbone");
                obj.bbConnectionPlatformPresent = 0;
                return
            end
                
            % Add a platform at the connection point with the backbone
            if obj.bbConnectionPlatformPresent
                if obj.H2Interconnect
                    obj.bbPlatform = Platform(grid, xIntrin,yIntrin, "H2toH2", "bb");
                else
                    if obj.scenario == "fullElectric"
                        obj.bbPlatform = Platform(grid, xIntrin,yIntrin, "EtoE", "bb");
                    else
                        obj.bbPlatform = Platform(grid, xIntrin,yIntrin, "EtoH2", "bb");
                    end
                end
            end
            
            % calculate distance from placed bbPlatform to two bbNodes
            pipe = chosenPipe;
            closestbbNodes = [];
            
            % find the corresponding bbNodes to the chosenPipe
            for i = 1:numel(grid.bbNodes)
                if any(grid.bbNodes(i).lines == pipe)
                    closestbbNodes(end+1) = i;
                end
            end
            
            % For each node, calculate the distance to shore
            numNodes = numel(closestbbNodes);
            min_length = nan;
            endNode = grid.shapes.connectionPoints(obj.onshoreConnection).bbNode;
            
            for i = 1:numNodes
                % Calculate distance from bblatfom 
                [coordX, coordY] = projfwd(grid.projection, grid.bbNodes(closestbbNodes(i)).coords(1), grid.bbNodes(closestbbNodes(i)).coords(2));
                dist = norm([coordX,coordY] - [end_coords(1), end_coords(2)]) / 1000;
                
                [lengthFromNode, nodesPassed, linesPassed] = grid.navigateBB(closestbbNodes(i), endNode);
                total_length = dist + lengthFromNode;
                    
                if isnan(min_length) || total_length < min_length
                    min_length = total_length;
                    obj.bbLength = min_length;
                    obj.bbLines = linesPassed;
                end
            end
            
            % Create a pipe/cable that's of the same length and properties of the
            % backbone
            if obj.scenario == "fullElectric"
                % Create cable
                transp = Cable(obj.hub2shoreCableVRating, obj.hub2shoreCableIRating, obj.hub2shoreCableFreq, obj.hub2shoreCableA, obj.hub2shoreCableCap);
            else
                % Create pipe
                transp = Pipe(obj.hub2shorePipeRadius, obj.hub2shorePipePressure, obj.hub2shorePipeEff, 50, 15);
            end
            
            transp.length = obj.bbLength;
            transp.bb = true;
            obj.bbTransport = transp;
        end
        
        function calculate_power(obj)
            % Rated power calculation
            for i = 1:numel(obj.turbines)
                obj.ratedPower = obj.ratedPower + obj.turbines(i).rating;
            end
            
            obj.inputPower = 0;
            % sum inputpowers of all cables/pipes. This is expected
            % averaged power coming from the turbines
            if obj.H2Interconnect
                lines = obj.pipes;
                outLines = obj.outPipes;
            else
                lines = obj.cables;
                if obj.scenario == "fullElectric" || obj.scenario == "electricToH2"
                    outLines = obj.outCables;
                else
                    outLines = obj.outPipes;
                end
            end
            
            numLines = numel(lines);
            
            for i = 1:numLines
                lines(i).calculate_power();
                obj.inputPower = obj.inputPower + lines(i).inputPower;
            end
            
            % Take care of any potential subplatform
            if obj.centralSubPresent
                % Take into account the efficieny of the subplatform
                subOutPower = obj.subPlatform.calculate_power(obj.inputPower);
                
                % Then the transportation losses via lines
                numOutlines = numel(outLines);
                OutLinePower = 0;
                
                for i = 1:numOutlines
                    if (subOutPower / numOutlines) > outLines(i).power_rating
                        error("Power exceeds cable rating")
                    end
                    outLines(i).set_power(subOutPower / numOutlines);
                    OutLinePower = OutLinePower + outLines(i).outputPower;
                end
            else
                numOutlines = numel(outLines);
                OutLinePower = 0;
                
                for i = 1:numOutlines
                    outLines(i).calculate_power();
                    OutLinePower = OutLinePower + outLines(i).outputPower;
                end
            end
            

            % Then the bbplatform losses
            if obj.bbConnectionPlatformPresent
                if obj.centralSubPresent
                    obj.bbPlatform.inPressure = obj.sub2mainPipePressure;
                    obj.bbPlatform.inVoltage = obj.sub2mainCableVRating;
                    obj.bbPlatform.outPressure = obj.hub2shorePipePressure;
                    obj.bbPlatform.outVoltage = obj.hub2shoreCableVRating;
                else
                    obj.bbPlatform.inPressure = obj.turb2subPipePressure;
                    obj.bbPlatform.inVoltage = obj.turb2subCableVRating;
                    obj.bbPlatform.outPressure = obj.hub2shorePipePressure;
                    obj.bbPlatform.outVoltage = obj.hub2shoreCableVRating;
                end
                bbInPower = obj.bbPlatform.calculate_power(OutLinePower);
                obj.outputPower = bbInPower;
            else
                obj.outputPower = OutLinePower;
            end
            
            % Then the lossen from bbPlatform to shore
            if obj.bbLength > 0
                obj.bbTransport.set_power(obj.outputPower);
                obj.shorePower = obj.bbTransport.outputPower;
            end
            
            obj.transport2shoreLosses = obj.inputPower - obj.shorePower;
            
            % Update properties
            if obj.H2Interconnect
                obj.pipes = lines;
                obj.outPipes = outLines;
            else
                obj.cables = lines;
                if obj.scenario == "fullElectric"
                    obj.outCables = outLines;
                else
                    obj.outPipes = outLines;
                end
            end
            
            obj.inputEnergy = obj.inputPower * 8760 / 1e9;
            obj.outputEnergy = obj.outputPower * 8760 / 1e9;
            obj.shoreEnergy = obj.shorePower * 8760 / 1e9;
            obj.capacityFactor = obj.inputPower / obj.ratedPower;
        end
        
        function cost = calculateCost(obj)
            global costsParams;
            global property;
            
            % Calculate total CAPEX
            % calculate turbine capex
            turbCapex = 0;
            for i = 1:numel(obj.turbines)
                % turb equipment + installation costs
                turbCapex = turbCapex + obj.turbines(i).calculateCost();
            end
            
            transpCapex = 0;
            % Calculate cables capex
            if numel(obj.cables) > 0
                for i = 1:numel(obj.cables)
                    transpCapex = transpCapex + obj.cables(i).calculateCost();
                end
            end
            if numel(obj.outCables) > 0
                for i = 1:numel(obj.outCables)
                    transpCapex = transpCapex + obj.outCables(i).calculateCost();
                end
            end
            % calculate pipes capex
            if numel(obj.pipes) > 0
                for i = 1:numel(obj.pipes)
                    transpCapex = transpCapex + obj.pipes(i).calculateCost();
                end
            end
            if numel(obj.outPipes) > 0
                for i = 1:numel(obj.outPipes)
                    transpCapex = transpCapex + obj.outPipes(i).calculateCost();
                end
            end
            % backbone capex
            if numel(obj.bbTransport) > 0
                for i = 1:numel(obj.bbTransport)
                    transpCapex = transpCapex + obj.bbTransport(i).calculateCost();
                end
            end
            
            % platform capex & opex
            platfCapex = 0;
            platfOpex = 0;
            if numel(obj.bbPlatform) > 0
                for i = 1:numel(obj.bbPlatform)
                    platfCapex = platfCapex + obj.bbPlatform(i).calculateCost();
                    platfOpex = platfOpex + obj.bbPlatform(i).OPEX;
                end
            end
            if numel(obj.subPlatform) > 0
                for i = 1:numel(obj.subPlatform)
                    platfCapex = platfCapex + obj.subPlatform(i).calculateCost();
                    platfOpex = platfOpex + obj.subPlatform(i).OPEX;
                end
            end
            
            total_equipm_install_cost = turbCapex + transpCapex + platfCapex;
            obj.CAPEXlist.turbines = turbCapex;
            obj.CAPEXlist.transport = transpCapex;
            obj.CAPEXlist.platforms = platfCapex;
            
            % calculate opex costs per component
            obj.OPEXlist.turbines = obj.CAPEXlist.turbines * costsParams.maintenance;
            obj.OPEXlist.transport = obj.CAPEXlist.transport * costsParams.maintenance;
            obj.OPEXlist.platforms = obj.CAPEXlist.platforms * costsParams.maintenance + platfOpex;

            % Calculate other costs
            managementCost = costsParams.management * (obj.ratedPower / 1e9);
            insuranceCost = total_equipm_install_cost * costsParams.insurance;
            
            obj.CAPEXlist.management = managementCost;
            obj.CAPEXlist.insurance = insuranceCost;
            
            obj.CAPEX = total_equipm_install_cost + managementCost + insuranceCost;
            
            % calculate LCOE
            obj.calculateLCOE();
        end
        
        function calculateLCOE(obj)
            global costsParams;
            
            energyProduction = obj.shoreEnergy * 1e3; % MWh
            
            % Calculate LCOE contribution of the turbines
            if numel(obj.turbines) > 0
                annuity_factor = (costsParams.WACC * (1 + costsParams.WACC)^obj.turbines(1).lifetime) / ((1+costsParams.WACC)^obj.turbines(1).lifetime - 1);
                obj.LCOElist.turbines = (obj.CAPEXlist.turbines * annuity_factor + obj.OPEXlist.turbines) * 1e6 / energyProduction;
            else
                obj.LCOElist.turbines = 0;
            end
            
            % Calculate LCOE contribution of the transport components
            if numel(obj.bbTransport) > 0
                annuity_factor = (costsParams.WACC * (1 + costsParams.WACC)^obj.bbTransport(1).lifetime) / ((1+costsParams.WACC)^obj.bbTransport(1).lifetime - 1);
                obj.LCOElist.transport = (obj.CAPEXlist.transport * annuity_factor + obj.OPEXlist.transport) * 1e6 / energyProduction;
            else
                obj.LCOElist.transport = 0;
            end
            
            % Calculate LCOE contribution of the platform components
            if numel(obj.bbPlatform) > 0
                annuity_factor = (costsParams.WACC * (1 + costsParams.WACC)^obj.bbPlatform.lifetime) / ((1+costsParams.WACC)^obj.bbPlatform.lifetime - 1);
                obj.LCOElist.platforms = (obj.CAPEXlist.platforms * annuity_factor + obj.OPEXlist.platforms) * 1e6 / energyProduction;
            else
                obj.LCOElist.platforms = 0;
            end
            
            % Calculate LCOE contribution of other costs using annuity
            % factor of the turbines
            annuity_factor = (costsParams.WACC * (1 + costsParams.WACC)^obj.turbines(1).lifetime) / ((1+costsParams.WACC)^obj.turbines(1).lifetime - 1);
            obj.LCOElist.other = ((obj.CAPEXlist.management + obj.CAPEXlist.insurance) * annuity_factor) * 1e6 / energyProduction;
            
            % Sum all individual LCOE components
            obj.LCOE = obj.LCOElist.turbines + obj.LCOElist.transport + obj.LCOElist.platforms + obj.LCOElist.other;
        end
        
        function cost = calculateCostOnLoc(obj)
            global costsParams;
            global property;
            % Calculate CAPEX without tranport to shore
            % calculate turbine capex
            turbCapex = 0;
            for i = 1:numel(obj.turbines)
                % turb equipment + installation costs
                turbCapex = turbCapex + obj.turbines(i).calculateCost();
            end
            
            transpCapex = 0;
            % Calculate inter-array cables capex
            if numel(obj.cables) > 0
                for i = 1:numel(obj.cables)
                    transpCapex = transpCapex + obj.cables(i).calculateCost();
                end
            end
            % calculate inter-array pipes capex
            if numel(obj.pipes) > 0
                for i = 1:numel(obj.pipes)
                    transpCapex = transpCapex + obj.pipes(i).calculateCost();
                end
            end
            
            total_equipm_install_cost = turbCapex + transpCapex;
            
            % Calculate other costs
            managementCost = costsParams.management * (obj.ratedPower / 1e9);
            insuranceCost = total_equipm_install_cost * costsParams.insurance;
            
            obj.CAPEX = total_equipm_install_cost + managementCost + insuranceCost;
            obj.CAPEXlist.turbines = turbCapex;
            obj.CAPEXlist.transport = transpCapex;
            obj.CAPEXlist.platforms = 0;
            obj.CAPEXlist.management = managementCost;
            obj.CAPEXlist.insurance = insuranceCost;
            
            % calculate opex
            obj.OPEX = costsParams.maintenance * total_equipm_install_cost;
            obj.OPEXlist.turbines = obj.CAPEXlist.turbines * costsParams.maintenance;
            obj.OPEXlist.transport = obj.CAPEXlist.transport * costsParams.maintenance;
            obj.OPEXlist.platforms = obj.CAPEXlist.platforms * costsParams.maintenance;
            
            % calculate LCOE/LCOH
            obj.calculateLCOE();
            
%             if obj.scenario == "H2inTurb"
%                 obj.H2Production = (energyProduction * 1e6) / property.H2specEnergy; % kg
%                 obj.LCOH = (obj.CAPEX * costsParams.annuityFactor * 1.03 + obj.OPEX) * 1e6 / obj.H2Production;
%             end
        end
        
        function mean = calc_mean_wind(obj)
            mean = 0;
            num = 0;
            for i = 1:numel(obj.turbines)
                mean = mean + obj.turbines(i).wind;
                num = num + 1;
            end
            mean = mean/num;
            obj.mean_wind = mean;
        end
        
        function mean = calc_mean_depth(obj)
            mean = 0;
            num = 0;
            for i = 1:numel(obj.turbines)
                mean = mean + obj.turbines(i).depth;
                num = num + 1;
            end
            mean = mean/num;
            obj.mean_depth = mean;
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
            for i = 1:numel(obj.outCables)
                if any(obj.outCables(i).xIntrin)
                    lenCable = max(size(obj.outCables(i).xIntrin));
                    xCoord = [];
                    yCoord = [];
                    for k = 1:lenCable
                        xCoord(end+1) = grid.X(obj.outCables(i).xIntrin(k), obj.outCables(i).yIntrin(k));
                        yCoord(end+1) = grid.Y(obj.outCables(i).xIntrin(k), obj.outCables(i).yIntrin(k));
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
            
            % Plot output pipes
            for i = 1:numel(obj.outPipes)
                if any(obj.outPipes(i).xIntrin)
                    lenPipe = max(size(obj.outPipes(i).xIntrin));
                    xCoord = [];
                    yCoord = [];
                    for k = 1:lenPipe
                        xCoord(end+1) = grid.X(obj.outPipes(i).xIntrin(k), obj.outPipes(i).yIntrin(k));
                        yCoord(end+1) = grid.Y(obj.outPipes(i).xIntrin(k), obj.outPipes(i).yIntrin(k));
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

