props;

% resolution is 7 * turbine diameter
desired_spacing = property.turbRadius * 14 / 1000;

if (exist('g') == 1)
    % g exists
    if (g.resolution == desired_spacing)
        disp("Grid already has the desired spacing.");
        return
    end
end

g = grid(desired_spacing);
% 
% g.add_tiff_data(property.maps.windSpeed, "wind", NaN);
% g.add_tiff_data(property.maps.windVariance, "windVar", NaN);
g.add_wind_data('GIS_data/mean_windspeed.mat');
g.add_tiff_data(property.maps.routeDensity, "traffic", -9999);
g.add_tiff_data(property.maps.bathymetry, "bath", 'bath');

g.add_nature_zones(property.shapes.natura2000, property.shapes.CDDA ); % has to be this order...
g.add_fishing_zones(property.shapes.fishing, property.fishing_threshold);
g.add_windfarm_zones(property.shapes.windfarms);
g.add_pipelines(property.shapes.pipelines);
g.add_landareas();
g.createBackboneGraph();
g.add_connectionPoints(property.shapes.onshoreConnections);

g.createGraph();

g.resetMask()