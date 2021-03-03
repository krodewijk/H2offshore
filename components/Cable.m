classdef Cable < Transport
    % Cable
    
    properties
        power_rating; % W
        voltage_rating; % V
        current_rating = 715; % A
        frequency;
        crossA; % mm2, cross section of cable
        capacitance; % uF / km, capacitance of cable
        
        resistivity; % Ohm * mm2 / m
        Rkm;
        
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
            
            % Calculate DC resistance / km
            obj.Rkm = obj.resistivity / obj.crossA * 1000;
            
            % Calculate power rating
            if obj.frequency > 0
                obj.power_rating = sqrt(3) * obj.current_rating * obj.voltage_rating * obj.power_factor;
            else
                obj.power_rating = obj.current_rating * obj.voltage_rating;
            end
                
            obj.power_capacity = obj.power_rating;
            
            obj.connected.turbines = Turbine.empty;
        end
        
        function calculate_power(obj)
            % Calculate input power (from turbines)
            numTurbs = numel(obj.connected.turbines);
            powerSum = 0;
            for i = 1:numTurbs
                powerSum = powerSum + obj.connected.turbines(i).power;
            end
            obj.inputPower = powerSum;
           
            % Calculate losses
            act_I = powerSum / obj.voltage_rating;
            
            % ohmic
            R = obj.resistivity * obj.length / 1000;
            Pohm = act_I ^2 * R;
            
            % dielectric
            C = obj.capacitance * obj.length / 1000; % uF
            Pdi = 2 * pi * obj.frequency * C * 10^-6 * obj.voltage_rating^2;
            
            % total lost power
            Ptot = Pohm + Pdi;
            
            obj.outputPower = obj.inputPower - Ptot;
        end
        
    end
end

