% convergence testing of sch_1d_cn()

% Case 1
idtype = 0;
vtype  = 0;
idpar  = [3];
vpar   = 0;
tmax   = 0.25;
lambda = 0.1;
lmin   = 6;
lmax   = 9;

[dpsi6rms, t6] = ctest_1d_func(tmax, 6, lambda, idtype, idpar, vtype, vpar);
[dpsi7rms, t7] = ctest_1d_func(tmax, 7, lambda, idtype, idpar, vtype, vpar);
[dpsi8rms, t8] = ctest_1d_func(tmax, 8, lambda, idtype, idpar, vtype, vpar);

figure
hold on
plot(t6, dpsi6rms,     'r-.')
plot(t7, 4*dpsi7rms,   'g-.')
plot(t8, 4*4*dpsi8rms, 'b-.')

xlabel("time")
ylabel("||d\psi^l||_2")
title("Time evolution of ||d\psi^l||_2, scaled by \rho = [1,4,4^2]. Case 1.")

legend("l=6", "l=7", "l=8")
hold off

[err6, t6] = ctest_1d_func_err(tmax, 6, lambda, idtype, idpar, vtype, vpar);
[err7, t7] = ctest_1d_func_err(tmax, 7, lambda, idtype, idpar, vtype, vpar);
[err8, t8] = ctest_1d_func_err(tmax, 8, lambda, idtype, idpar, vtype, vpar);

figure
hold on
plot(t6, dpsi6rms,     'r-.')
plot(t7, 4*dpsi7rms,   'g-.')
plot(t8, 4*4*dpsi8rms, 'b-.')

xlabel("time")
ylabel("||E(\psi^l)||_2")
title("Time evolution of ||E(\psi^l)||_2, scaled by \rho = [1,4,4^2]. Case 1.")

legend("l=6", "l=7", "l=8")
hold off

% ----------------------------------------

% Case 2
idtype = 1;
vtype  = 0;
idpar  = [0.50 0.075 0.0];
vpar   = 0;
tmax   = 0.01;
lambda = 0.01;
lmin   = 6;
lmax   = 9;

[dpsi6rms, t6] = ctest_1d_func(tmax, 6, lambda, idtype, idpar, vtype, vpar);
[dpsi7rms, t7] = ctest_1d_func(tmax, 7, lambda, idtype, idpar, vtype, vpar);
[dpsi8rms, t8] = ctest_1d_func(tmax, 8, lambda, idtype, idpar, vtype, vpar);

figure
hold on
plot(t6, dpsi6rms,     'r-.')
plot(t7, 4*dpsi7rms,   'g-.')
plot(t8, 4*4*dpsi8rms, 'b-.')

xlabel("time")
ylabel("||d\psi^l||_2")
title("Time evolution of ||d\psi^l||_2, scaled by \rho = [1,4,4^2]. Case 2.")

legend("l=6", "l=7", "l=8")
hold off
