function add_wind_data(obj, fileName)
    wind_lat = ncread(fileName, 'lat');
    wind_lon = ncread(fileName, 'lon');
    windspeeds = ncread(fileName, 'wspeed');

    % Take the correct height and time
    winds = windspeeds(:,:,1)-100;

    % Calculate limits of the map
    wind_latlim = double([min(wind_lat(:)), max(wind_lat(:))]);
    wind_lonlim = double([min(wind_lon(:)), max(wind_lon(:))]);

    % Convert the wind map to the already defined grid
    obj.data.wind = obj.convert_to_grid(wind_lat, wind_lon, winds);
end

function add_pipelines(obj, fileName)
    fprintf(1,' Adding pipelines\n');
    S = shaperead(fileName, 'UseGeoCoords', true);
    dim = size(S);
    dim = max(dim);

    gridded = nan(size(obj.lat));

    fprintf(1,' \nProgress:    ');
    % Add intrinsic coordinates to shape and add to grid
    for i = 1:dim % for each line
        innerDim = numel(S(i).Lat);
        pipeSize = str2double(S(i).SIZE_IN);

        % Create nodes field
        S(i).nodes = [];

        for k = 1:innerDim % for each coordinate in line make sure the nodes are within the grid

            % Add first point to grid, otherwise fit it first
            cur_lat = S(i).Lat(k);
            cur_lon = S(i).Lon(k);

            if isnan(cur_lat) || isnan(cur_lon)
                continue;
            end

            % Check if current coord is the last one
            if k ~= innerDim
                next_lat = S(i).Lat(k+1);
                next_lon = S(i).Lon(k+1);
            end

            if  k == innerDim || isnan(next_lat) || isnan(next_lon)
                next_lat = S(i).Lat(k-1);
                next_lon = S(i).Lon(k-1);
            end

            % Check if current coord is outside the grid
            if cur_lat < obj.latlim(1) || cur_lat > obj.latlim(2) || cur_lon < obj.lonlim(1) || cur_lon > obj.lonlim(2)
                % cur coordinate is outside grid

                % Calculate gradient
                grad = (next_lat - cur_lat) / (next_lon - cur_lon);

                % Calculate intersections with limits
                dlon = [(obj.latlim(1) - cur_lat); (obj.latlim(2) - cur_lat)] / grad;
                dlat = [(obj.lonlim(1) - cur_lon); (obj.lonlim(2) - cur_lon)] * grad;

                intersects = [obj.latlim(1), dlon(1) + cur_lon; obj.latlim(2), dlon(2) + cur_lon];
                intersects = [intersects; cur_lat + dlat(1), obj.lonlim(1); cur_lat + dlat(2), obj.lonlim(2)];

                point_dist = nan;
                % Check if intersec coords is between limits of map and
                % limits of pipeline
                for m = 1:4
                    cur_intersect = intersects(m,:);
                    if cur_intersect(1) >= obj.latlim(1) && cur_intersect(1) <= obj.latlim(2) && cur_intersect(2) >= obj.lonlim(1) && cur_intersect(2) <= obj.lonlim(2)
                        % Current intersect point lies in grid
                        if cur_intersect(1) >= min(cur_lat, next_lat) && cur_intersect(1) <= max(cur_lat, next_lat) && cur_intersect(2) >= min(cur_lon, next_lon) && cur_intersect(2) <= max(cur_lon, next_lon)
                            % Also in between min and max of line
                            % This is the new point
                            new_dist = norm([cur_lat, cur_lon] - [cur_intersect(1), cur_intersect(2)]);
                            if new_dist < point_dist || isnan(point_dist)
                                point_dist = new_dist;
                                new_lat = cur_intersect(1);
                                new_lon = cur_intersect(2);
                            end
                        end
                    end
                end

                cur_lat = new_lat;
                cur_lon = new_lon;
            end

            % Add the nodes to the grid graph
            [xIntrin, yIntrin] = obj.geo2intrin(cur_lat, cur_lon);

            if isnan(xIntrin) || isnan(yIntrin)
                % Still outside scope
                continue;
            end
            cur_node = obj.intrin2node(xIntrin, yIntrin);

            if k > 1
                prev_node = S(i).nodes(end);
                path = shortestpath(obj.graph, prev_node, cur_node);
                S(i).nodes = [S(i).nodes, path];
            else
                S(i).nodes(end+1) = cur_node;
            end

        %disp(k)
        end
        fprintf(1,'\b\b\b\b%3.0f%%',(i/dim*100));
        %disp(i);
    end

    fprintf(1,' \n');
    obj.shapes.pipelines = S;
end

function new_map = convert_to_grid(obj, lat_data, lon_data, data)
    % "interpolate" new map by assigning defined
    % coordinates to nearest new data
    [rows, cols] = size(obj.lat);
    new_map = nan([rows, cols]); % initial')ize new map with nans

    % Calculate limits of the new map
    data_latlim = double([min(lat_data(:)), max(lat_data(:))]);
    data_lonlim = double([min(lon_data(:)), max(lon_data(:))]);

    % Iterate over defined coordinates
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
            [min_distance, col_ind] = min(min_row);
            min_coords = [row_ind(col_ind), col_ind];

            new_map(x, y) = data(min_coords(1),min_coords(2));
        end
        disp([num2str(x), ' / ', num2str(rows), ' | ', num2str(x/rows*100), ' %'])
    end
end

function coordMap = make_coordMap(obj, lat_data, lon_data)
    disp('Creating coordinate map');
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

function createGraph
            [rows, cols] = size(obj.lat);
            numNodes = rows * cols;

            curNode = 1;

            A = sparse(numNodes, numNodes);

            fprintf(1,' \nProgress:    ');
            for y = 1:rows
                for x = 1:cols
                    % left node
                    if (x > 1)
                        leftNode = curNode - 1;
                        dist = norm([obj.X(y, x - 1), obj.Y(y, x - 1)] - [obj.X(y,x), obj.Y(y,x)]);
                        A(curNode, leftNode) = dist;
                    end
                    % right node
                    if (x < cols)
                        rightNode = curNode + 1;
                        dist = norm([obj.X(y, x + 1), obj.Y(y, x + 1)] - [obj.X(y,x), obj.Y(y,x)]);
                        A(curNode, rightNode) = dist;
                    end
                    % bottom node
                    if (y > 1)
                        bottomNode = curNode - cols;
                        dist = norm([obj.X(y - 1, x), obj.Y(y - 1, x)] - [obj.X(y,x), obj.Y(y,x)]);
                        A(curNode, bottomNode) = dist;
                    end
                    % top node
                    if (y < rows)
                        topNode = curNode + cols;
                        dist = norm([obj.X(y + 1, x), obj.Y(y + 1, x)] - [obj.X(y,x), obj.Y(y,x)]);
                        A(curNode, topNode) = dist;
                    end
                    curNode = curNode + 1;
                end
                fprintf(1,'\b\b\b\b%3.0f%%',(y/rows*100));
            end
end