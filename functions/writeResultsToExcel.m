function writeResultsToExcel(farm)
    data = [];

    % intrinsic coordinates
    descr(1) = "Matlab x coordinate";
    data(1) = farm.start_x;

    descr(end+1) = "Matlab y coordinate";
    data(end+1) = farm.start_y;

    % mean wind
    descr(end+1) = "Mean wind speed, 10m height (m/s)";
    data(end+1) = farm.mean_wind;

    descr(end+1) = "Mean wind speed, hub height (m/s)";
    data(end+1) = farm.turbines(1).translate_wind(farm.mean_wind);

    % mean depth
    descr(end+1) = "Mean depth";
    data(end+1) = abs(farm.mean_depth);

    % # turbines
    descr(end+1) = "Number of turbines";
    data(end+1) = numel(farm.turbines);

    % Length of backbone
    descr(end+1) = "Length of backbone to shore (km)";
    data(end+1) = farm.bbLength;

    % Pressure of backbone
    descr(end+1) = "Pressure at end of backbone (bar)";
    if farm.scenario == "H2inTurb"
        data(end+1) = farm.bbTransport.outPressure;
    else
        data(end+1) = 0;
    end

    % Full load hours
    descr(end+1) = "Full-load hours";
    data(end+1) = farm.capacityFactor;

    % Rated power
    descr(end+1) = "Farm rated power (MW)";
    data(end+1) = farm.ratedPower / 1e6; % MW

    % Yield turbine bottom
    descr(end+1) = "Yield (GWh/year)";
    data(end+1) = farm.inputEnergy; % GWh / year

    % Yield feeded in backbone
    descr(end+1) = "Yield supplied to backbone (GWh/year)";
    data(end+1) = farm.outputEnergy; % GWh / year

    % Yield arriving at shore
    descr(end+1) = "Yield at shore (GWh/year)";
    data(end+1) = farm.shoreEnergy; % GWh / year

    % Compression energy
    descr(end+1) = "Compression energy needed (MWh/year)";
    data(end+1) = farm.comprEnergy; % MWh

    % Mean turbine cost
    descr(end+1) = "Mean turbine cost (excluding transportation to shore) (M€)";
    data(end+1) = mean(arrayfun(@(x) x.costOnLoc, farm.turbines(:)));

    % Mean cost per turbine (including transportation)
    descr(end+1) = "Mean cost per turbine (including transportation to shore) (M€)";
    data(end+1) = mean(arrayfun(@(x) x.costOffshore, farm.turbines(:)));

    % Total cost of windfarm excluding transportation
    descr(end+1) = "Total cost of windfarm excluding transportation (M€)";
    data(end+1) = farm.costOnLoc;

    % Total cost of windfarm including transportation
    descr(end+1) = "Total cost of windfarm including transportation to shore (M€)";
    data(end+1) = farm.costOnshore;

    % Compression costs
    descr(end+1) = "Compression costs (M€)";
    data(end+1) = farm.comprPrice;

    % LCOE
    descr(end+1) = "LCOE on location (€/MWh)";
    data(end+1) = farm.LCOEOnLoc; % EUR / MWh
    descr(end+1) = "LCOE on shore (€/MWh)";
    data(end+1) = farm.LCOEOnshore; % EUR / MWh

    % LCOH
    descr(end+1) = "LCOE on location (€/kg)";
    data(end+1) = farm.LCOHOnLoc; % EUR / kg
    descr(end+1) = "LCOE on shore (€/kg)";
    data(end+1) = farm.LCOHOnshore; % EUR / kg

    str_tab = table(descr');
    writetable(str_tab, "output/matlab_results.xlsm", 'Range', 'B4','WriteVariableNames',false);

    if farm.scenario == "fullElectric"
        tab = table(data(1:end-2)');
        writetable(tab,"output/matlab_results.xlsm", 'Range', 'C4','WriteVariableNames',false);
    else
        tab = table(data');
        writetable(tab,'output/matlab_results.xlsm', 'Range', 'D4','WriteVariableNames',false);
    end

end