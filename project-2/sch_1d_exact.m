function [psi_exact] = sch_1d_exact(x, t, m)
% compute exact solution to 1D schrodinger equation with initial condition
%  psi(x,0) = sin(m*pi*x)
nx = length(x);
nt = length(t);

x_array = repmat(x,nt,1);
t_array = repmat(transpose(t),1,nx);

psi_exact = exp(-1i*m^2*pi^2*t_array).*sin(m*pi*x_array);
end

