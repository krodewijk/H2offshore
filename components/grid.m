classdef grid < handle
    %GRID Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        latlim = [50.5, 61.5]; % vertical angle
        lonlim = [-3.2, 8.2]; % horizontal angle
        projection = projcrs(23095); % ED50 North Sea projection
        %projection = projcrs(28992);
        resolution = 2.5; % km
        
        % geo coordinates
        lat;
        lon;
        
        % projected cartesian coordinates in meters
        X; 
        Y;
        
        % graph system
        graph;
        
        cellReference;
        data = struct;
        shapes = struct;
        mask;
        
        turbines;
    end
    
    methods
        function obj = grid(resolution)
            % Set object properties
            global property;
            obj.latlim = property.latlim;
            obj.lonlim = property.lonlim;     
            obj.resolution = resolution;
            obj.projection = projcrs(property.projSystem);
            
            %latlim_distance = ((obj.latlim(2) - obj.latlim(1)) / 360 * 2 * pi) * earthRadius('km');
            %lonlim_distance = ((obj.lonlim(2) - obj.lonlim(1)) / 360 * 2 * pi) * earthRadius('km');
            
            % Calculate km between lat limits
            [coordX1, coordY1] = projfwd(obj.projection, obj.latlim(1), obj.lonlim(1));
            [coordX2, coordY2] = projfwd(obj.projection, obj.latlim(2), obj.lonlim(1));
            
            latlim_distance = norm([coordX1, coordY1] - [coordX2, coordY2]) / 1000;
            
            % Calculate km betwoon lon limits
            [coordX1, coordY1] = projfwd(obj.projection, obj.latlim(1), obj.lonlim(1));
            [coordX2, coordY2] = projfwd(obj.projection, obj.latlim(1), obj.lonlim(2));
            
            lonlim_distance = norm([coordX1, coordY1] - [coordX2, coordY2]) / 1000;
            
            latres = round(latlim_distance / obj.resolution);
            lonres = round(lonlim_distance / obj.resolution);
            
            [obj.lat, obj.lon] = meshgrat(obj.latlim, obj.lonlim, [latres, lonres]);
            
            obj.cellReference = georefcells(obj.latlim, obj.lonlim, size(obj.lon));
            
            % Create projected XY grid
            [obj.X, obj.Y] = projfwd(obj.projection, obj.lat, obj.lon);
            
            % Create mask with ones
            obj.mask = true(size(obj.lat));
        end
        
        function add_fishing_zones(obj, fileName, min_MW)
            fprintf(1,' Adding Fishing zones\n');
            S = shaperead(fileName, 'Attributes',{'mw_fshn'}, 'UseGeoCoords', true, 'Selector',...
           {@(v1) (v1 >= min_MW),'mw_fshn'});
       
            obj.shapes.fishing = S;
            
            fishing_grid = zeros(size(obj.lat));
            
            fprintf(1,' Progress:    ');
            for i = 1:max(size(S))
                fishing_grid = or(fishing_grid, inpolygon(obj.lat, obj.lon, S(i).Lat, S(i).Lon));
                fprintf(1,'\b\b\b\b%3.0f%%',(i/max(size(S))*100));
            end
            fprintf(1,'\n');
            obj.data.fishing = fishing_grid;
        end
        
        function add_windfarm_zones(obj, fileName)
            fprintf(1,' Adding windfarm zones\n');
            S = shaperead(fileName,'UseGeoCoords', true);
       
            obj.shapes.windfarms = S;
            
            windfarm_grid = zeros(size(obj.lat));
            
            fprintf(1,' Progress:    ');
            for i = 1:max(size(S))
                windfarm_grid = or(windfarm_grid, inpolygon(obj.lat, obj.lon, S(i).Lat, S(i).Lon));
                fprintf(1,'\b\b\b\b%3.0f%%',(i/max(size(S))*100));
            end
            fprintf(1,'\n');
            obj.data.windfarm = windfarm_grid;
        end
        
        function add_landareas(obj)
            fprintf(1,' Adding landareas\n');
            S = shaperead('landareas','UseGeoCoords', true);
       
            obj.shapes.landareas = S;
            
            landareas_grid = zeros(size(obj.lat));
            
            fprintf(1,' Progress:    ');
            for i = 1:max(size(S))
                landareas_grid = or(landareas_grid, inpolygon(obj.lat, obj.lon, S(i).Lat, S(i).Lon));
                fprintf(1,'\b\b\b\b%3.0f%%',(i/max(size(S))*100));
            end
            fprintf(1,'\n');
            obj.data.landareasMask = landareas_grid;
        end
        
        function add_nature_zones(obj, fileName, fileName2)
            fprintf(1,' Adding Nature zones\n');
            S1 = shaperead(fileName, 'UseGeoCoords', true);
            S2 = shaperead(fileName2, 'UseGeoCoords', true);
            
            % Two bad apples in data...
            S1(101) = [];
            S1(101) = [];
       
            obj.shapes.nature1 = S1;
            obj.shapes.nature2 = S2;
            
            nature_grid = zeros(size(obj.lat));
            
            fprintf(1,' Progress Shape 1:    ');
            for i = 1:max(size(S1))
                nature_grid = or(nature_grid, inpolygon(obj.lon, obj.lat, S1(i).Lon, S1(i).Lat));
                fprintf(1,'\b\b\b\b%3.0f%%',(i/max(size(S1))*100));
                %disp(i)
            end
            
            fprintf(1,' \nProgress Shape 2:    ');
            for i = 1:max(size(S2))
                nature_grid = or(nature_grid, inpolygon(obj.lat, obj.lon, S2(i).Lat, S2(i).Lon));
                fprintf(1,'\b\b\b\b%3.0f%%',(i/max(size(S2))*100));
            end
            
            fprintf(1,'\n');
            obj.data.nature = nature_grid;
        end
        
        function add_tiff_data(obj, fileName, dataName, convNan)
            % convNan traffic => -9999
            
            fprintf(1,' Adding tif data %s\n', dataName);
            [data, R] = readgeoraster(fileName);
            
            % Convert the data map to the already defined grid and add to
            % properties
            gridded = geointerp(data, R, obj.lat, obj.lon);
            
            if(convNan == 'bath')
                gridded(gridded==max(gridded)) = NaN;
            else
                gridded(gridded==convNan) = NaN;
            end
            
            obj.data.(dataName) = gridded;
        end
        
        function add_wind_data(obj, fileName)
            disp('Adding wind data.');
            load(fileName, 'mean_wind');
            
            % remove coordinates of data that fall outside of defined
            % grid to speed up the progress
            disp('Remove data outside of defined grid');
            [newRows, newCols] = size(mean_wind.lat);
            for x = 1:newRows
                for y = 1:newCols
                    if (mean_wind.lat(x,y) < obj.latlim(1)) || (mean_wind.lat(x,y) > obj.latlim(2)) || (mean_wind.lon(x,y) < obj.lonlim(1)) || (mean_wind.lon(x,y) > obj.lonlim(2))
                        % Lies outside limits
                        % Remove data
                        mean_wind.lat(x,y) = NaN;
                        mean_wind.lon(x,y) = NaN;
                        mean_wind.speed(x,y) = NaN;
                        mean_wind.var(x,y) = NaN;
                        mean_wind.hist(x,y,:) = NaN;
                    end
                end
            end
            
            % CUT AWAY ALL THE NANs (sorry for the code...)
            foundData = false;
            
            disp('Remove empty rows');
            removeRows = [];
            removeCols = [];
            % Remove nan rows
            x = 1;
            while x < numel(mean_wind.lat(:,1))
                if min(isnan(mean_wind.lat(x,:)))
                    % Whole row consists of nans
                    % Remove row
                    removeRows(end+1) = x;
                end
                x = x+1;
            end
            
            mean_wind.lat(removeRows,:) = [];
            mean_wind.lon(removeRows,:) = [];
            mean_wind.speed(removeRows,:) = [];
            mean_wind.var(removeRows,:) = [];
            mean_wind.hist(removeRows,:,:) = [];
            
            disp('Remove empty columns');
            % Remove nan columns
            y = 1;
            while y < numel(mean_wind.lat(1,:))
                if min(isnan(mean_wind.lat(:,y)))
                    % Whole column consists of nans
                    % Remove column
                    removeCols(end+1) = y;
                end
                y = y+1;
            end
            
            mean_wind.lat(:,removeCols) = [];
            mean_wind.lon(:,removeCols) = [];
            mean_wind.speed(:,removeCols) = [];
            mean_wind.var(:,removeCols) = [];
            mean_wind.hist(:,removeCols,:) = [];
            
            % Create coordinate mapping to grid object
            coordMap = obj.make_coordMap(mean_wind.lat, mean_wind.lon);
            
            [rows, cols] = size(obj.lat);
            
            % define windDist
            obj.data.windDist = nan(rows, cols, numel(mean_wind.bins));
            % Add data to grid
            fprintf(1,' \nAdding data to object:    ');
            for x = 1:rows
                for y = 1:cols
                    data_coords = squeeze(coordMap(x,y,:));
                    if isnan(data_coords(1)) || isnan(data_coords(2))
                        obj.data.wind(x,y) = nan;
                        obj.data.windVar(x,y) = nan;
                        obj.data.windDist(x,y,:) = nan(100,1);
                    else
                        obj.data.wind(x,y) = mean_wind.speed(data_coords(1), data_coords(2));
                        obj.data.windVar(x,y) = mean_wind.var(data_coords(1), data_coords(2));
                        obj.data.windDist(x,y,:) = mean_wind.hist(data_coords(1), data_coords(2),:);
                    end 
                end
                fprintf(1,'\b\b\b\b%3.0f%%',(x/rows*100));
            end
            % Also add bins of wind distribution
            obj.data.windDistBins = mean_wind.bins';
            
            clearvars mean_wind;
        end
        
        function coordMap = make_coordMap(obj, lat_data, lon_data)
            disp('Redefining coordinate map (this will take a while...)');
            % "interpolate" new map by assigning defined
            % coordinates to nearest new data
            [rows, cols] = size(obj.lat);
            
            coordMap = nan([rows, cols, 2]); % initialize coordinate map --> each intrinsic coordinate of coordMap is the same as grid object and has the corresponding new intrin coords of the new data

            % Calculate limits of the new map
            data_latlim = double([min(lat_data(:)), max(lat_data(:))]);
            data_lonlim = double([min(lon_data(:)), max(lon_data(:))]);

            % Iterate over defined coordinates
            fprintf(1,' \nProgress:    ');
            for x = 1:rows
                for y = 1:cols
                    % Update current lat & lon
                    cur_lat = obj.lat(x,y);
                    cur_lon = obj.lon(x,y);

                    % Make sure the current coordinates are within the new
                    % data map limits
                    if (cur_lat < data_latlim(1)) || (cur_lat > data_latlim(2))
                        continue
                    elseif (cur_lon < data_lonlim(1)) || (cur_lon > data_lonlim(2))
                        continue
                    end

                    % Calculate distances from current point to defined
                    % grid
                    distances = distance(cur_lat, cur_lon, lat_data, lon_data);

                    % find indices of minimum value in distances matrix
                    [min_row, row_ind] = min(distances); % returns minimum row
                    [~, col_ind] = min(min_row);
                    min_coords = [row_ind(col_ind), col_ind];

                    coordMap(x, y, :) = min_coords;
                end
                fprintf(1,'\b\b\b\b%3.0f%%',(x/rows*100));
            end
        end
        
        function add_connectionPoints(obj, fileName)
            fprintf(1,' Adding onshore connection points\n');
            S = shaperead(fileName,'UseGeoCoords', true);
            
            % Add intrinsic coordinates + node to struct
            for i = 1:numel(S)
                [xIntrin, yIntrin] = obj.geo2intrin(S(i).Lat, S(i).Lon);
                S(i).xIntrin = xIntrin;
                S(i).yIntrin = yIntrin;
                node = obj.intrin2node(xIntrin, yIntrin);
                S(i).node = node;
            end
            
            obj.shapes.connectionPoints = S;
        end

        function createGraph(obj)
            fprintf(1,' \nCreating graph in grid.');
            [rows, cols] = size(obj.lat);
            numNodes = rows * cols;

            curNode = 1;

            A = sparse(numNodes, numNodes);

            fprintf(1,' \nProgress:    ');
            for x = 1:rows
                for y = 1:cols
                    % curNode
                    % sweep from left to right, up and down around 3 sized square around middle 
                    for x_in = (x-1):(x+1)
                        for y_in = (y-1):(y+1)
                            % check if x and y are within bound
                            if(x_in > 0 && y_in > 0 && x_in <= rows && y_in <= cols)
                                endNode = obj.intrin2node(x_in, y_in);
                                dist = norm([obj.X(x_in, y_in), obj.Y(x_in, y_in)] - [obj.X(x,y), obj.Y(x,y)]);
                                A(curNode, endNode) = dist;
                            end
                        end
                    end
                    curNode = curNode + 1;
                end
                fprintf(1,'\b\b\b\b%3.0f%%',(x/rows*100));
            end

            obj.graph = graph(A);
        end
        
        function add_pipelines(obj, fileName)
            fprintf(1,' Adding pipelines\n');
            S = shaperead(fileName, 'UseGeoCoords', true);
            dim = size(S);
            dim = max(dim);
            
            fprintf(1,' \nProgress:    ');
            % Add intrinsic coordinates to shape and add to grid
            for i = 1:dim % for each line
                innerDim = numel(S(i).Lat);
                pipeSize = str2double(S(i).SIZE_IN);
                
                % Create nodes field
                S(i).nodes = [];
                
                for k = 1:innerDim % for each coordinate in line
                    
                    % Add first point to grid, otherwise fit it first
                    cur_lat = S(i).Lat(k);
                    cur_lon = S(i).Lon(k);
                    
                    if isnan(cur_lat) || isnan(cur_lon)
                        continue;
                    end
                        
                    % Add the nodes to the grid graph
                    [xIntrin, yIntrin] = obj.geo2intrin(cur_lat, cur_lon);

                    if isnan(xIntrin) || isnan(yIntrin)
                        % Outside scope
                        continue;
                    end
                    cur_node = obj.intrin2node(xIntrin, yIntrin);
                    
                    S(i).nodes(end+1) = cur_node;

                %disp(k)
                end
                fprintf(1,'\b\b\b\b%3.0f%%',(i/dim*100));
                %disp(i);
            end
            
            fprintf(1,' \n');
            obj.shapes.pipelines = S;
        end
        
        function num = intrin2node(obj, x, y)
            dim = size(obj.lat);
            num = dim(2) * (x-1) + y;
        end

        function [x,y] = node2intrin(obj, node)
            dim = size(obj.lat);
            x = ceil(node / dim(2));
            y = mod(node, dim(2));
        end
        
        function [xIntrin,yIntrin] = geo2intrin(obj, lat, lon)
            % Make sure the current coordinates are within the
            % new map limits
            if (lat < obj.latlim(1)) || (lat > obj.latlim(2) || isnan(lat))
                xIntrin = NaN;
                yIntrin = NaN;
            elseif (lon < obj.lonlim(1)) || (lon > obj.lonlim(2) || isnan(lon))
                xIntrin = NaN;
                yIntrin = NaN;
            else
                % Is within bounds
                xIntrin = round(obj.cellReference.latitudeToIntrinsicY(lat));
                yIntrin = round(obj.cellReference.longitudeToIntrinsicX(lon));
            end
        end
        
        function [lat, lon] = intrin2geo(obj, xIntrin, yIntrin)
            lon = obj.cellReference.intrinsicXToLongitude(yIntrin);
            lat = obj.cellReference.intrinsicYToLatitude(xIntrin);
        end
        
        function path = shortestRoute(obj, start_x, start_y, end_x, end_y)
            path = shortestpath(obj.graph, obj.intrin2node(start_x, start_y), obj.intrin2node(end_x, end_y));
        end
        
        function resetMask(obj)
            global property;
            obj.mask = true(size(obj.lat));
            %obj.mask = obj.mask & ~obj.data.nature & ~obj.data.fishing & (obj.data.traffic < property.traffic_threshold);

            % Add land areas
            obj.mask = obj.mask & ~obj.data.landareasMask;
        end
        
        function plot_pipes(obj)
            coastLines = load("coastlines");
            cities = shaperead('worldcities', 'UseGeoCoords', true);
            citiesLat = zeros(max(size(cities),1));
            citiesLon = zeros(max(size(cities)),1);
            
            % Convert shape format to useable array
            for i = 1:max(size(cities))
                citiesLat(i) = cities(i).Lat;
                citiesLon(i) = cities(i).Lon;
            end
            
            % Plot pipes2shore
            numPipes = numel(obj.shapes.pipelines);
            for i = 1:numPipes
                lenCable = numel(obj.shapes.pipelines(i).xIntrin);
                xCoord = [];
                yCoord = [];
                
                [xCoord, yCoord] = projfwd(obj.projection, obj.shapes.pipelines(i).Lat, obj.shapes.pipelines(i).Lon);

%                 for k = 1:lenCable
%                     if ~isnan(obj.shapes.pipelines(i).xIntrin(k)) && ~isnan(obj.shapes.pipelines(i).yIntrin(k))
%                         xCoord(end+1) = obj.X(obj.shapes.pipelines(i).xIntrin(k), obj.shapes.pipelines(i).yIntrin(k));
%                         yCoord(end+1) = obj.Y(obj.shapes.pipelines(i).xIntrin(k), obj.shapes.pipelines(i).yIntrin(k));
%                 
%                     end
%                 end
                plot(xCoord, yCoord,'b', 'LineWidth',5);
                hold on;
            end

            % Plot surroundings (coastlines + cities)
            [xProjCoast, yProjCoast] = projfwd(obj.projection, coastLines.coastlat, coastLines.coastlon);
            [xProjCity, yProjCity] = projfwd(obj.projection, citiesLat, citiesLon);
            plot(xProjCoast, yProjCoast);
            plot(xProjCity, yProjCity, 'o', 'MarkerSize', 10, 'LineWidth', 5);
            
            
            xlim([-2e5, 7.5e5]);
            ylim([5.4e6, 7.4e6]);
        end
        
        function plot(obj,dataName)
            mapaxis = axesm('lambertstd', 'MapLatLimit', obj.latlim, 'MapLonLimit', obj.lonlim);
            mapshow(obj.X, obj.Y, double(obj.data.(dataName)), 'DisplayType', 'surface');
            
            colormap;
            caxis([0, 15]);

        end
    end
end

