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
plat.inputPower = 1e9;
plat.inPressure = 100;
plat.outPressure = 300;
plat.H2toH2PowerTransfer();
compr_out = plat.compPower;
compr_energy = compr_out * 3600