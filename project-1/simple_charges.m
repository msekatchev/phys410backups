nc = 13;
tmax = 1000;
level = 15;
gamma = 1;
epsec = 1.0e-5;

r0 = random_r0(nc);

[t, r, v, v_ec] = charges(r0, tmax, level, gamma, epsec);