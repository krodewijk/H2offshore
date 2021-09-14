classdef Pipe < Transport
    % Pipe class for H2 transportation
    
    properties
        radius; % m
        inPressure; % bar input pressure
        outPressure; % bar output pressure
        pipe_efficiency;
        inTemp; % input temperature degC
        outTemp; % output temperature degC
        maxFlow; % m3 / day
        actualFlow; % m3 / day
        energy_loss; % in % / km
        
        power_rating; % in W
        
        H2density; % kg/m3 @ 15degC, 1 bar
        H2specEnergy; % Wh / kg
        H2compressibility; 
        H2gravity;
        H2normalGasTemp; % K
        baseTemp; % K
        basePress; % bar
        maxV;
    end
    
    methods
        function obj = Pipe(radius, pressure, pipe_eff, inTemp, outTemp)
            global property;
            global costsParams;
            
            obj.H2density = property.H2density;
            obj.H2specEnergy = property.H2specEnergy;
            obj.H2compressibility = property.H2compressibility; 
            obj.H2normalGasTemp = property.H2normalGasTemp; % K
            obj.baseTemp = property.baseTemp; % K
            obj.basePress = property.basePress; % bar
            obj.H2gravity = property.H2gravity;
            obj.maxV = property.maxV;
            obj.bbBasePerc = property.BBbaseload;
            obj.lifetime = costsParams.cableLifetime;
            
            obj.radius = radius;
            obj.inPressure = pressure;
            obj.pipe_efficiency = pipe_eff;
            obj.inTemp = inTemp;
            obj.outTemp = outTemp;
            
            obj.maxFlow = pi * obj.radius^2 * obj.maxV * 3600 * 24; % m3/day standard
            
            % Calculate power rating
            actual_density = (obj.inPressure * 0.75) * obj.H2density; % estimation 3/4 the input pressure is taken as mean pipe pressure
            
            obj.power_rating = actual_density * obj.maxFlow * obj.H2specEnergy * 1/24;  
            %obj.power_rating = obj.H2density * obj.maxFlow * obj.H2specEnergy * 1/24;  %at base conditions
            obj.power_capacity = obj.power_rating;
            obj.bbBasePower = obj.bbBasePerc * obj.power_rating;
            
            obj.connected.turbines = Turbine.empty;
        end
        
        function calculate_power(obj)
            numTurbs = numel(obj.connected.turbines);
            % only calculate input pwer if turbines are connected,
            % otherwise set the inputPower porperty before calling this
            % method
            if numTurbs > 0
                powerSum = 0;
                for i = 1:numTurbs
                    powerSum = powerSum + obj.connected.turbines(i).power;
                end
                obj.inputPower = powerSum;
            end
            
            % Calculate pressure drop
            % Calculate actual gas flow
            % obj.inputPower = 120e6;
            
            if obj.bb
                % calculate gasflow of base power + inputPower
                gasFlowTot = (obj.inputPower+obj.bbBasePower) * (obj.basePress * obj.H2density)^-1 * obj.H2specEnergy^-1 * 24;
                gasFlowBase = obj.bbBasePower * (obj.basePress * obj.H2density)^-1 * obj.H2specEnergy^-1 * 24;
                
                obj.outPressure = obj.pressure_drop(gasFlowTot) + (obj.inPressure - obj.pressure_drop(gasFlowBase));
            else
                gasFlow = obj.inputPower * (obj.basePress * obj.H2density)^-1 * obj.H2specEnergy^-1 * 24;
                obj.actualFlow = gasFlow;
            
                obj.outPressure = obj.pressure_drop(gasFlow);
            end
            
            if obj.outPressure < 0
                error("Pressure drop is too big");
            end
                
            %obj.outputPower = obj.outPressure / obj.inPressure * obj.inputPower;
            obj.outputPower = obj.inputPower;
            %obj.energy_loss = (obj.outputPower / obj.inputPower) / obj.length * 100;
            obj.energy_loss = 0;
        end
        
        function outPressure = pressure_drop(obj, gasFlow)
            % Calculate pressure drop
            outP2 = (obj.inPressure * 100)^2 - obj.H2gravity * obj.H2normalGasTemp * obj.length * obj.H2compressibility * (267.13 * gasFlow * obj.pipe_efficiency^(-1) * ((obj.basePress*100)/obj.baseTemp) * (2*obj.radius*1000)^(-2.667))^2;
            
            outPressure = sqrt(outP2) / 100; % in bar
        end
        
        function cost = calculateCost(obj)
            global costsParams;
            if obj.bb
                cost = costsParams.hub2shorePipe * obj.length;
            else
                cost = costsParams.turb2subPipe * obj.length;
            end
            obj.CAPEX = cost;
        end
    end
end

