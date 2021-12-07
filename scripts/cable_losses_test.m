%% CABLE
% A1 --> 775
c = Cable(525, 5625, 0, 3000, 0);

c.inputPower = 1e9;
c.length = 100;
c.calculate_power();

c_loss = c.inputPower - c.outputPower;
c_perc = c_loss / c.inputPower * 100;

%% PIPE
% r1 --> 0.1524 m
% r2 -->0.6095 m
p = Pipe(0.6095, 80, 0.95, 15,15);
p.inputPower = 1e9:1e9:13e9;
p.length = 1000;
p.calculate_power();
press_loss = p.inPressure - p.outPressure;
press_perc = press_loss / p.inPressure * 100;

%% Compressor
plat = Platform(g, 200, 200, "H2toH2", "bb");
plat.inputPower = 1e9; % W
H2_energy = plat.inputPower * 8760; % Wh/year
H2_production = H2_energy / property.H2specEnergy; %kg

N = 200;
end_pressure = 800;
compr_out = zeros(N,1);
compr_ratio = zeros(N,1);
pressures = zeros(N,1);
for i = 1:N
    pressures(i) = end_pressure / N * i;
    plat.inPressure = 1;
    plat.outPressure = pressures(i);
    plat.H2toH2PowerTransfer();
    compr_out(i) = plat.compPower * 8760 / 1000; % kWh/year
    compr_ratio(i) = compr_out(i) / H2_production; % kWh / kg
end

plot(pressures, compr_ratio, 'k');
title("Energy required for hydrogen compression");
ylabel("kWh/kg");
xlabel("Pressure (bar)");