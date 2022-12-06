function [psi_exact] = sch_2d_exact(x, y, t, mx, my)
% compute exact solution to 2D schrodinger equation with initial condition
%  psi(x,0) = sin(mx*pi*x)sin(my*pi*y)
nx = length(x);
ny = length(y);
nt = length(t);

psi_exact = zeros(nt,nx,ny, 'like', 1j);

x_array = repmat(x,nt,1);
y_array = repmat(y,nt,1);

for ti = [1:nt]
    psi_exact(ti,:,:) = exp(-1i*(mx^2+my^2)*pi^2*t(ti)).*sin(mx*pi*x).*sin(my*pi*y.');
end
end

