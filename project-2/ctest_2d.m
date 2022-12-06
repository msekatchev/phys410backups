% convergence testing of sch_2d_adi()

% Case 1
idtype = 0;
vtype  = 0;
idpar  = [2, 3];
vpar   = 0;
tmax   = 0.05;
lambda = 0.05;
lmin   = 6;
lmax   = 9;

disp("6")
[dpsi6rms, t6] = ctest_2d_func(tmax, 6, lambda, idtype, idpar, vtype, vpar);
disp("7")
[dpsi7rms, t7] = ctest_2d_func(tmax, 7, lambda, idtype, idpar, vtype, vpar);
disp("8")
[dpsi8rms, t8] = ctest_2d_func(tmax, 8, lambda, idtype, idpar, vtype, vpar);

figure
hold on
plot(t6, dpsi6rms,     'r-.')
plot(t7, 4*dpsi7rms,   'g-.')
plot(t8, 4*4*dpsi8rms, 'b-.')

xlabel("time")
ylabel("||d\psi^l||_2")
title("Evolution of ||d\psi^l||_2, scaled by \rho = [1,4,4^2]. Case 1.")

legend("l=6", "l=7", "l=8")
hold off




[err6, t6] = ctest_2d_func_err(tmax, 6, lambda, idtype, idpar, vtype, vpar);
[err7, t7] = ctest_2d_func_err(tmax, 7, lambda, idtype, idpar, vtype, vpar);
[err8, t8] = ctest_2d_func_err(tmax, 8, lambda, idtype, idpar, vtype, vpar);

figure
hold on
plot(t6, dpsi6rms,     'r-.')
plot(t7, 4*dpsi7rms,   'g-.')
plot(t8, 4*4*dpsi8rms, 'b-.')

xlabel("time")
ylabel("||E(\psi^l)||_2")
title("Evolution of ||E(\psi^l)||_2, scaled by \rho = [1,4,4^2]. Case 1.")

legend("l=6", "l=7", "l=8")
hold off



