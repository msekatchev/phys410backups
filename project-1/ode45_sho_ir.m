% ode45_sho_ir: Integrates equations of motion for simple harmonic
% oscillator using ODE45 and evaluates independent residuals using
% a finite difference approximation
% Integrate on the domain 0 <= t <= 3 pi, with initial conditions
%
% y_1(0) = x(0) = 0
% y_2(1) = v(0) = 1
%
% corresponding to the exact solution
%
% y(t) = sin(t)
% Integrate with default parameters, except use an uniform mesh of
% output times of length 2^10 + 1 and a stringent error tolerance
tol = 1.0e-10;
options = odeset('AbsTol', [tol, tol], 'RelTol', tol);
[tout yout] = ode45(@fcn_sho, linspace(0.0,3.0*pi,1025), [0.0 1.0], ...
options);
% Compute and plot three levels of scaled independent residuals
figure(1); clf; figure(2); clf
yir = yout(:,1);
tir = tout;
for ll = 1 : 3
dt = tir(2) - tir(1);
ir = zeros(1,length(yir));
for j = 2 : length(yir) - 1
ir(j) = (yir(j-1) - 2.0 * yir(j) + yir(j+1)) / dt^2 + yir(j);
end
% Plot unscaled residuals
figure(1);
hold on;
plot(tir, ir);
% Plot scaled residuals
figure(2);
hold on;
plot(tir, 4.0^(3 - ll) * ir);
yir = yir(1:2:end);
tir = tir(1:2:end);
end;