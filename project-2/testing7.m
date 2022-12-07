tmax = 0.05;
level = 6;
lambda = 0.05;
idtype = 0;
idpar = [2,2];
vtype = 0;
vpar = 0;

[x y t psi psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);


[Y,X] = meshgrid(y,x);

animate_thing_2df(1, X, Y, psire, 1);
