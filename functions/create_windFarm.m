% Search in 10x10 grid and create wind farms
wfs = Windfarm.empty;

[rows, cols] = size(g.lat);
wf_index = 1;

fprintf(1,'\nProgress:    ');
for x = 1:10:(rows-1)
    for y = 1:10:(cols-1)
        wfs(wf_index) = Windfarm(g, x, y, wf_index);
        wf_index = wf_index + 1;
        %disp([num2str(x), ', ', num2str(y)]);
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(x/rows*100));
end
fprintf(1,'\n');