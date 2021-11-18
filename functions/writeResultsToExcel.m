function writeResultsToExcel(farm)
    data = [];

    % intrinsic coordinates
    data(1) = farm.start_x;
    data(2) = farm.start_y;

    % mean wind
    data(3) = farm.mean_wind;
    data(4) = farm.turbines(1).translate_wind(farm.mean_wind);

    % mean depth
    data(5) = farm.mean_depth;

    % # turbines
    data(6) = numel(farm.turbines);

    % Length of backbone
    data(7) = farm.bbLength;

    % Full load hours
    data(8) = farm.capacityFactor;

    % Rated power
    data(9) = farm.ratedPower / 1e6; % MW

    % Yield turbine bottom
    data(10) = farm.inputEnergy; % GWh / year

    % Yield feeded in backbone
    data(11) = farm.outputEnergy; % GWh / year

    % Yield arriving at shore
    data(12) = farm.shoreEnergy; % GWh / year

    % LCOE
    data(13) = farm.LCOEOnLoc; % EUR / MWh
    data(14) = farm.LCOEOnshore; % EUR / MWh

    % LCOH
    data(15) = farm.LCOHOnLoc; % EUR / kg
    data(16) = farm.LCOHOnshore; % EUR / kg

    if farm.scenario == "fullElectric"
        tab = table(data(1:end-2)');
        writetable(tab,"output/matlab_results_electric.csv");
    else
        tab = table(data');
        writetable(tab,'output/matlab_results_h2.csv');
    end

end