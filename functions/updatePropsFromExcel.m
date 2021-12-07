function [property, costTables] = updatePropsFromExcel(fileName)
% Load the parameters from the excel file into the property struct
    props;

    excel = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'D7:E49');
    criticalRows = [2;3;4;7;8;9;10;11;14;15;19;16;17;18;20;21;22;23;24;25;41;29;28;30;31;32;33;36;37;40;41];
    
    if any(isnan(excel.Var2(criticalRows)))
        % a critical value in excel is missing
        % give an error
        error("A critical parameter in excel is missing. Please check the excel sheet.");
    end
    % READ AND STORE ALL VALUES IN EXCEL TO CORRESPONDING PROPERTIES
    % Simulation scenario
    property.scenario = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'K4:K5').Variables;
    property.scenario = property.scenario{1};

    if ~(property.scenario == "fullElectric" || property.scenario == "H2inTurb")
        error("Ïncorrect simulation scenario in Excel file")
    end
    % Farm size and power
    property.farmDim = [excel.Var2(2), excel.Var2(3)];
    property.farmPower = excel.Var2(4); % GW

    % Turbine parameters
    property.turbRating = excel.Var2(7) * 1e6; % W
    property.turbHubHeight = excel.Var2(8); % m
    property.turbRadius = excel.Var2(9) / 2; % m
    property.turbMaxWind = excel.Var2(10); % m/s
    property.turbMinWind = excel.Var2(11); % m/s

    % Cable parameters
    property.turb2subCableVRating = excel.Var2(14); % kV
    property.turb2subCableIRating = excel.Var2(15); % A (current rating)
    property.turb2subCableFreq = excel.Var2(19); % Hz
    property.turb2subCableA = excel.Var2(16); % mm2 (cross-section)
    property.turb2subCableCap = excel.Var2(17); % uF / km (capacitance)
    property.turb2subCableRho = excel.Var2(18); % Ohm * mm2/m (copper resistivity)

    property.hub2shoreCableVRating = excel.Var2(20); % kV (voltage rating)
    property.hub2shoreCableIRating = excel.Var2(21); % A (current rating)
    property.hub2shoreCableFreq = excel.Var2(25); % if cable is AC or DC (true = AC, false = DC)
    property.hub2shoreCableA = excel.Var2(22); % mm2 (cross-section)
    property.hub2shoreCableCap = excel.Var2(23); % uF / km (capacitance)
    property.hub2shoreCableRho = excel.Var2(24); % Ohm * mm2/m (copper resistivity)

    % Pipeline parameters
    property.turb2subPipePressure = excel.Var2(41); % bar (pressure rating)
    property.turb2subPipeRadius = excel.Var2(29) * 0.0254 / 2; % m (pipe radius)
    property.turb2subPipeEff = excel.Var2(28); % pipeline efficiency

    property.main2bbPipePressure = excel.Var2(32); % bar (pressure rating)
    property.main2bbPipeRadius = excel.Var2(31) * 0.0254 / 2; % m (pipe radius)
    property.main2bbPipeEff = excel.Var2(30); % pipeline efficiency
    property.BBbaseload = excel.Var2(33); % 60% baseload

    % H2 compression parameters
    property.platform.comprEff = excel.Var2(36); % percent
    property.platform.comprAdiaEff = excel.Var2(37); % percent
    property.platform.comprOutPress = property.main2bbPipePressure; % bar

    % Electrolyser parameters
    property.platform.electrolyserEff = excel.Var2(40); % percent
    property.platform.electrolyserOutPress = excel.Var2(41); % Bar


    %% READ AND STORE ALL DEPTH AND DISTANCE DEPENDENT COST PARAMETERS FROM EXCEL
    warning('off','all');
    costTables.fixed.grid = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F57:W76', 'ReadRowNames', true, 'ReadVariableNames', true);
    costTables.fixed.turbine = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F78:W97', 'ReadRowNames', true, 'ReadVariableNames', true);
    costTables.fixed.h2toshore = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F99:W118', 'ReadRowNames', true, 'ReadVariableNames', true);
    costTables.fixed.h2onloc = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F141:W160', 'ReadRowNames', true, 'ReadVariableNames', true);
    
    costTables.float.grid = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F162:W181', 'ReadRowNames', true, 'ReadVariableNames', true);
    costTables.float.turbine = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F183:W202', 'ReadRowNames', true, 'ReadVariableNames', true);
    costTables.float.h2toshore = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F204:W223', 'ReadRowNames', true, 'ReadVariableNames', true);
    costTables.float.h2onloc = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F246:W265', 'ReadRowNames', true, 'ReadVariableNames', true);
    
    costTables.breakevens = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'G268:W268');
    costTables.distances = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'G57:W57');
    costTables.depths = readtable(fileName, 'Sheet', 'Matlab Export', 'Range', 'F58:F76');
    warning('on','all')
end