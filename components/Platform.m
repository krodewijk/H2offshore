classdef Platform < handle
    % Central Hub of every cluster. Receives power from windfarms and
    % output power to shore
    % Platform has electrolyser+compressor/ACDC Converter 
    
    properties
        % Location of platform in grid coordinates
        xIntrin;
        yIntrin;
        node;
        
        kind; % EtoH2, EtoE, H2toH2
        floating;
        depth;
        
        inputPower = 0;
        outputPower = 0;
        
        electrolyserEff;
        converterEff;
        comprEff;
    end
    
    methods
        function obj = Platform(grid, xIntrin, yIntrin)
            global property;
            
            % Constructor
            obj.xIntrin = xIntrin;
            obj.yIntrin = yIntrin;
            obj.node = grid.intrin2node(xIntrin, yIntrin);
            
            obj.electrolyserEff = property.platform.electrolyserEff / 100;
            obj.converterEff = property.platform.transfEff / 100;
            obj.comprEff = property.platform.comprEff / 100;
            
            obj.floating = property.platform.floating;
            obj.depth = grid.data.bath(xIntrin, yIntrin);
            
            scenario = property.scenario;
            
            if scenario == "H2inTurb"
                obj.kind = "H2toH2";
            elseif scenario == "fullElectric"
                obj.kind = "EtoE";
            elseif scenario == "electricToH2"
                obj.kind = "EtoH2";
            else
                error("Scenario not set correctly.")
            end
            
            % Flip mask 
            grid.mask(xIntrin, yIntrin) = false;
        end
        
        function outPower = calculate_power(obj, inputPower) 
            obj.inputPower = inputPower;
            if obj.kind == "H2toH2"
                outPower = obj.H2toH2PowerTransfer(inputPower);
            elseif obj.kind == "EtoE"
                outPower = obj.EtoEPowerTransfer(inputPower);
            elseif obj.kind == "EtoH2"
                outPower = obj.EtoH2PowerTransfer(inputPower);
            else
                error("Scenario not set correctly.")
            end
            
            obj.outputPower = outPower;
        end
        
        function outPower = H2toH2PowerTransfer(obj, inputPower)
            % Only compressor needed for H2 to H2 transport
            outPower = obj.comprEff * inputPower;
        end
        
        function outPower = EtoEPowerTransfer(obj, inputPower)
            % Only converter needed
            outPower = obj.converterEff * inputPower;
        end
        
        function outPower = EtoH2PowerTransfer(obj, inputPower)
            % Electrolyser + compressor needed
            outPower = obj.converterEff * obj.electrolyserEff * inputPower;
        end
    end
end

