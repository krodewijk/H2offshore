dataSize = size(mean_wind.speed);
rows = dataSize(1);
cols = dataSize(2);
structIndex = 1;

for i = 1:rows
    for k = 1:cols
        DOWA(structIndex).Geometry = 'Point';
        DOWA(structIndex).Lon = double(mean_wind.lon(i,k));
        DOWA(structIndex).Lat = double(mean_wind.lat(i,k));
        DOWA(structIndex).speed = double(mean_wind.speed(i,k));
        structIndex = structIndex + 1;
    end
end

shapewrite(DOWA, "DOWA.shp");