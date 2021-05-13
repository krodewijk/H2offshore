classdef Platform < handle
    % Any platform that may be needed
    % output power to shore
    % Platform has electrolyser/compressor/ACDC Converter 
    
    properties
        % Location of platform in grid coordinates
        xIntrin;
        yIntrin;
        node;
        
        kind; % EtoH2, EtoE, H2toH2
        location; % substation ("sub") or backbone compressor station ("bb")
        floating;
        depth;
        
        inputPower = 0;
        outputPower = 0;
        
        electrolyserEff;
        converterEff;
        comprEff;
        comprAdiaEff;
        
        % H2toH2 parameters
        inPressure;
        outPressure;
        H2compressibility;
        H2specHeatRatio;
        gasTemp;
        H2density;
        H2specEnergy;
        basePress;
        
        %EtoE paramters
        inVoltage;
        outVoltage;
        
        CAPEX; % M-EUR
    end
    
    methods
        function obj = Platform(grid, xIntrin, yIntrin, kind, loc)
            global property;
            
            % Constructor
            obj.xIntrin = xIntrin;
            obj.yIntrin = yIntrin;
            obj.node = grid.intrin2node(xIntrin, yIntrin);
            
            obj.electrolyserEff = property.platform.electrolyserEff / 100;
            obj.converterEff = property.platform.transfEff / 100;
            obj.comprEff = property.platform.comprEff / 100;
            obj.comprAdiaEff = property.platform.comprAdiaEff / 100;
            obj.H2compressibility = property.H2compressibility;
            obj.H2specHeatRatio = property.specHeatRatio;
            obj.gasTemp = property.baseTemp;
            obj.basePress = property.basePress;
            obj.H2density = property.H2density;
            obj.H2specEnergy = property.H2specEnergy;
            
            obj.floating = property.platform.floating;
            obj.depth = grid.data.bath(xIntrin, yIntrin);
            obj.location = loc;
            
%             if scenario == "H2inTurb"
%                 obj.kind = "H2toH2";
%             elseif scenario == "fullElectric"
%                 obj.kind = "EtoE";
%             elseif scenario == "electricToH2"
%                 obj.kind = "EtoH2";
%             else
%                 error("Scenario not set correctly.")
%             end
            obj.kind = kind;
            
            if obj.location == "sub"
                % H2toH2
                obj.inPressure = property.turb2subPipePressure;
                obj.outPressure = property.sub2mainPipePressure;
                
                % EtoE
                obj.inVoltage = property.turb2subCableVRating;
                obj.outVoltage = property.sub2hubCableVRating;
                
            elseif obj.location == "bb"
                % H2toH2
                obj.inPressure = property.sub2mainPipePressure;
                obj.outPressure = property.main2bbPipePressure;
                
                % EtoE
                obj.inVoltage = property.sub2hubCableVRating;
                obj.outVoltage = property.hub2shoreCableVRating;
            end
            
            % Flip mask 
            grid.mask(xIntrin, yIntrin) = false;
        end
        
        function outPower = calculate_power(obj, inputPower) 
            obj.inputPower = inputPower;
            if obj.kind == "H2toH2"
                outPower = obj.H2toH2PowerTransfer();
            elseif obj.kind == "EtoE"
                outPower = obj.EtoEPowerTransfer();
            elseif obj.kind == "EtoH2"
                outPower = obj.EtoH2PowerTransfer();
            else
                error("Scenario not set correctly.")
            end
            
            obj.outputPower = outPower;
        end
        
        function outPower = H2toH2PowerTransfer(obj)
            % Only compressor needed for H2 to H2 transport
            flowRate = obj.inputPower * (obj.basePress * obj.H2density)^-1 * obj.H2specEnergy^-1 * 24 / 1e6; % Mm3 / day

            compPower = 4063.9 * (obj.H2specHeatRatio / (obj.H2specHeatRatio-1)) * flowRate * obj.gasTemp * obj.H2compressibility * (1/obj.comprAdiaEff) * ((obj.outPressure / obj.inPressure)^((obj.H2specHeatRatio - 1)/obj.H2specHeatRatio)-1);
            
            outPower = obj.inputPower - compPower;
        end
        
        function outPower = EtoEPowerTransfer(obj)
            % Only converter needed
            outPower = obj.converterEff * obj.inputPower;
        end
        
        function outPower = EtoH2PowerTransfer(obj)
            % Electrolyser + compressor needed
            outPower = obj.electrolyserEff * obj.inputPower;
        end
        
        function cost = calculateCost(obj)
            global costsParams;
            
            platformCost = costsParams.platformEquipm; % Million
            installCost = costsParams.platformInstallation; % Million
            
            % Component cost
            % Compressor/transformer??
            % Electrolyser
            if obj.kind == "EtoH2"
                compCost = costsParams.H2ProductionCAPEX * obj.inputPower / 1e9; % M-EUR
            else
                compCost = 0;
            end
            
            cost = platformCost + installCost + compCost;
            
            obj.CAPEX = cost;
        end
    end
end

