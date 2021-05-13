classdef Cable < Transport
    % Cable
    
    properties
        power_rating; % W
        voltage_rating; % V
        current_rating = 715; % A
        frequency;
        crossA; % mm2, cross section of cable
        capacitance; % uF / km, capacitance of cable
        dielectricLossFactor;
        
        resistivity; % Ohm * mm2 / m
        Rkm;
        act_I;
        
        power_factor = 0.95;
    end
    
    methods
        function obj = Cable(kv_rating, i_rating, freq, A, cap)
            global property;
            
            obj.voltage_rating = kv_rating * 1e3;
            obj.frequency = freq;
            obj.current_rating = i_rating;
            obj.crossA = A;
            obj.capacitance = cap;
            obj.resistivity = property.turb2subCableRho;
            obj.dielectricLossFactor = property.dielectricLossFactor;
            obj.bbBasePerc = property.BBbaseload;
            
            % Calculate DC resistance / km
            obj.Rkm = obj.resistivity / obj.crossA * 1000;
            
            % Calculate power rating
            if obj.frequency > 0
                obj.power_rating = sqrt(3) * obj.current_rating * obj.voltage_rating * obj.power_factor;
            else
                obj.power_rating = obj.current_rating * obj.voltage_rating;
            end
                
            obj.power_capacity = obj.power_rating;
            obj.bbBasePower = obj.bbBasePerc * obj.power_rating;
            
            obj.connected.turbines = Turbine.empty;
        end
        
        function outP = calculate_power(obj)
            numTurbs = numel(obj.connected.turbines);
            if numTurbs > 0
                % Calculate input power (from turbines)
                powerSum = 0;
                for i = 1:numTurbs
                    powerSum = powerSum + obj.connected.turbines(i).power;
                end
                obj.inputPower = powerSum;
            end
            
            if obj.bb
                % Add baseload power/current to the cable
                Ploss = obj.cable_losses(obj.inputPower + obj.bbBasePower) - obj.cable_losses(obj.bbBasePower);
            else
                Ploss = obj.cable_losses(obj.inputPower);
            end
            
            obj.outputPower = obj.inputPower - Ploss;

            outP = obj.outputPower;
            
            if obj.outputPower < 0
                warning("Too much losses in cable!")
            end
        end
        
        function loss = cable_losses(obj, inPower)
            act_I = inPower / obj.voltage_rating;
            if obj.frequency > 0
                act_I = act_I / sqrt(3);
            end
            
            % ohmic
            R = obj.Rkm * obj.length;
            Pohm = act_I ^2 * R;
            
            % dielectric
            C = obj.capacitance * obj.length; % uF
            Pdi = 2 * pi * obj.frequency * C * 10^-6 * obj.voltage_rating^2 * obj.dielectricLossFactor;
            
            % total lost power
            loss = Pohm + Pdi;
        end
        
        function cost = calculateCost(obj)
            global costsParams;
            if obj.bb
                cableCost = costsParams.hub2shoreCosts * obj.length;
            else
                cableCost = costsParams.turb2subCosts * obj.length;
            end
            installCost = costsParams.cableInstallation * obj.length;
            cost = cableCost + installCost;
            obj.CAPEX = cost;
        end
    end
end

