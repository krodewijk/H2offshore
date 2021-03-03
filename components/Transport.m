classdef Transport < handle
    % Network class
    
    properties
        % Start and end coordinates in intrins coord
        startNode;
        endNode;
        
        % Intrinsic coordinates where pipe will go
        xIntrin;
        yIntrin;
        
        length=0; % km
        
        nodes;

        inputPower = 0;
        outputPower = 0;
        
        connected = struct;
        power_capacity; % W capacity left in transport
    end

    
    methods
        function success = add_turbine(obj, grid, Turbine)
            % Check if there is capactiy left on the cable for another
            % turbine
            if (obj.power_capacity > Turbine.rating)
                % Add turbine to list
                obj.connected.turbines(end+1) = Turbine;
                obj.power_capacity = obj.power_capacity - Turbine.rating;
                success = true;
            else
                success = false;
                return;
            end
            if(max(size(obj.connected.turbines)) > 1)
                obj.update_connections(grid);
            elseif (numel(obj.connected.turbines) == 1)
                obj.nodes = obj.connected.turbines.node;
            end
            obj.calculate_power();
            obj.calc_length(grid);
        end
        
        function remove_turbine(obj, grid, Turbine)
            % Check if turbine is in list
            if (min(size(obj.connected.turbines(obj.connected.turbines == Turbine))) > 0)
                % Turbine is in list
                obj.connected.turbines(obj.connected.turbines == Turbine) = [];
                obj.power_capacity = obj.power_capacity + Turbine.rating;
                % Update connections
                if(max(size(obj.connected.turbines)) > 1)
                    obj.update_connections(grid);
                end
            else
                error("Turbine was not connected");
            end
            obj.calculate_power();
        end
        
        function add_connection(obj, grid, newNode)
            % Check if a cable already exist
            if (max(size(obj.nodes)) > 0)
                % A path already exist
                start = obj.nodes(1);
                if numel(obj.nodes) > 1
                    eind = obj.nodes(2);
                else
                    eind = start;
                end
                
                % Check which distance is shortest
                start2new = distances(grid.graph, start, newNode);
                eind2new = distances(grid.graph, eind, newNode);
                if (start2new < eind2new)
                    % start is closest
                    newPath = shortestpath(grid.graph, start, newNode);
                    obj.nodes = [obj.nodes, newPath];
                    obj.nodes = unique(obj.nodes, 'stable');
                else
                    % eind is closest
                    newPath = shortestpath(grid.graph, eind, newNode);
                    obj.nodes = [obj.nodes, newPath];
                    obj.nodes = unique(obj.nodes, 'stable');
                end
            else
                % No path already exist
                obj.nodes = [newNode];
            end
            obj.update_intrinCoords(grid);
            obj.calculate_power();
            obj.calc_length(grid);
        end
        
        function len = calc_length(obj, grid)
            numNodes = numel(obj.nodes)-1;
            len = 0;
            
            for i = 1:numNodes
                dist = distances(grid.graph, obj.nodes(i), obj.nodes(i+1));
                len = len + dist;
            end
            
            obj.length = len;
        end
        
        function update_connections(obj, grid)
            % Sort Turbine array on nodeNumber
            numTurbs = max(size(obj.connected.turbines));
            sorted_turbs = Turbine.empty;
            turbs_path = Turbine.empty;
            obj.nodes = [];
            
            [~, ind] = sort([obj.connected.turbines.node]);
            
            sorted_turbs = obj.connected.turbines(ind);
            
            % First one in turbs path is first one of sorted list
            turbs_path(1) = sorted_turbs(1);
            added_turbs = [sorted_turbs(1).node]; % nodes of added turbs

            % Get sorted nodes
            i = 2;
            while (i <= numTurbs)
                % Check if next turbine is already added
                if any(added_turbs(:) == sorted_turbs(i).node)
                    i = i+1;
                    continue;
                end
                
                % Get xy coords
                [xCur, yCur] = grid.node2intrin(turbs_path(end).node);
                [xNext, yNext] = grid.node2intrin(sorted_turbs(i).node);
                % Check if next turbine is next to current (hor or vert)
                % and if not already added to path
                if ((abs(xCur - xNext) > 1) || (abs(yCur - yNext) > 1)) || ((abs(xCur - xNext) == 1) && (abs(yCur - yNext) == 1))
                    % Next turbine not directly next to current
                    % turbine
                    closestTurb = sorted_turbs(i);
                    closestDist = norm([closestTurb.xIntrin, closestTurb.yIntrin] - [xCur, yCur]);
                    % Find closest turbine
                    for k = i:numTurbs
                        [coordX, coordY] = grid.node2intrin(sorted_turbs(k).node);
                        dist = norm([coordX, coordY] - [xCur, yCur]);
                        if (dist < closestDist && ~any(added_turbs(:) == sorted_turbs(k).node))
                            closestTurb = sorted_turbs(k);
                            closestDist = dist;
                        end
                    end
                    % Add closest Turb to path
                    turbs_path = [turbs_path, closestTurb];
                    added_turbs = [added_turbs, closestTurb.node];
                    i = i - 1;
                else
                    % Next turbine is next to current
                    % Add next Turb to path
                    turbs_path = [turbs_path, sorted_turbs(i)];
                    added_turbs = [added_turbs, sorted_turbs(i).node];
                end
                
                % Get shortest path from previous turbine
                path = shortestpath(grid.graph, turbs_path(end-1).node, turbs_path(end).node);
                % Add path to cable nodes
                obj.nodes = [obj.nodes, path];
                i = i + 1;
            end
            
            % Remove duplicate nodes
            obj.nodes = unique(obj.nodes, 'stable');
            
            obj.update_intrinCoords(grid);
            
            obj.startNode = obj.nodes(1);
            obj.endNode = obj.nodes(end);
        end
        
        function update_intrinCoords(obj, grid)
            obj.xIntrin = [];
            obj.yIntrin = [];
            for i = 1:max(size(obj.nodes))
                % Add path to instrinsic coords
                [xIntrin, yIntrin] = grid.node2intrin(obj.nodes(i));
                obj.xIntrin(end+1) = xIntrin;
                obj.yIntrin(end+1) = yIntrin;
            end
        end
    end
end

