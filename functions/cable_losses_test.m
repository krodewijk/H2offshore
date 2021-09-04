c = Cable(525, 4952, 0, 3000, 0);

c.inputPower = 200e6;
c.length = 100;
c.calculate_power();

c_loss = (c.outputPower - c.inputPower) / c.inputPower;

p = Pipe(0.61, 80, 0.95, 0,0);
p.inputPower = 200e6;
p.length = 100;
p.calculate_power();
p_loss = (p.outputPower - p.inputPower) / p.inputPower;