clf
nc = 12;
tmax = 10;
level = 12;
gamma = 1;
epsec = 1.0e-5;

r0 = random_r0(nc);

[t, r, v, v_ec] = charges(r0, tmax, level, gamma, epsec);

plot(t,v, "-r")
xlabel("time")
ylabel("total potential energy V(t)")
title("Time evolution of potential for 12-charge calculation")
message = sprintf('tmax = 10\nlevel = 12\ngamma = 1')
annotation('textbox', [.2 .5 .2 .3], 'String', message,'FitBoxToText','on')
