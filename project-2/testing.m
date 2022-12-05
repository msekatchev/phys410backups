sch_1d_cn(0.5, 2, 0.25, 0, [0], 0, [0])


tmax = 0.25
level = 2
lambda = 0.1
idtype = 0
idpar = [0]
vtype = 0
vpar = 0

% initialize discretization mesh
nx = 2^level + 1;
x = linspace(0,1,nx);
dx = x(2) - x(1);
dt = lambda * dx;
nt = round(tmax / dt) + 1;

t   = [0: nt-1] * dt;
psi = zeros(nt,nx);
% initialize psi depending on IC type
if idtype == 0
    m = idpar(1);
    psi(1,:) = sin(m*pi*x);
else
    x0 = idpar(1);
    delta = idpar(2);
    p = idpar(3);
    psi(1,:) = exp(1i*p*x)*exp(  -((x-x0)/delta).^2  );
end

% initialize diagonals for tridiagonal matrix A
du = zeros(nx,1); % upper diagonal  "+"
dm = zeros(nx,1); % middle diagonal "0"
dl = zeros(nx,1); % lower diagonal  "-"
f  = zeros(nx,1); % RHS

du = 1/2 * 1/(dx)^2 * ones(nx,1);
dm = (1i/dt - 1/(dx)^2) * ones(nx,1);
dl = du;

% set up middle diagonal based on potential


du(2)    = 0;
dm(1)    = 1; % why?
dm(nx)   = 1; % why?
dl(nx-1) = 0;

% define sparse matrix:
A = spdiags([dl dm du], -1:1, nx, nx);

