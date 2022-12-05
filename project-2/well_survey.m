% use sch_1d_cn to compute the excess fractional probability for the 1D
% Schrodinger equation with a potential well.

tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.40, 0.075, 0.0];
vtype = 1;
xmin = 0.6;
xmax = 0.8;
x1 = 0.6;
x2 = 0.8;

num_points = 251;

Fe = zeros(1,num_points);
V0 = -exp(linspace(2,10,num_points));

for n = [1:num_points]
    disp(n)
    vpar = [xmin, xmax, V0(n)];
    Fe(n) = excess_frac_prob(tmax, level, lambda, idtype, idpar, vtype, vpar, x1, x2);
end

figure;

hold on;

plot(log(abs(V0)),log(Fe), "-");
xlabel('log$(|V_0|)$','Interpreter','latex')
ylabel('log$(\bar{F}_e(0.6,0.8))$','Interpreter','latex')
title("Well survey - excess fractional probability vs well depth")

hold off;
