function [x t psi psire psiim psimod prob v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Solves 1D time-dependent Schrodinger equation via the Crank-Nichelson method
% Inputs
%
% tmax: Maximum integration time
% level: Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parameters
%
% Outputs
%
% x: Vector of x coordinates [nx]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx]
% psire Array of computed psi_re values [nt x nx]
% psiim Array of computed psi_im values [nt x nx]
% psimod Array of computed sqrt(psi psi*) values [nt x nx]
% prob Array of computed running integral values [nt x nx]
% v Array of potential values [nx]

% initialize discretization mesh
nx = 2^level + 1;
x = linspace(0,1,nx);
dx = x(2) - x(1);
dt = lambda * dx;
nt = round(tmax / dt) + 1;

t   = [0: nt-1] * dt;
psi = zeros(nt,nx);

prob = zeros(nt,nx);

% initialize psi depending on IC type
if idtype == 0
    m = idpar(1);
    psi(1,:) = sin(m*pi*x);
else
    x0 = idpar(1);
    delta = idpar(2);
    p = idpar(3);
    psi(1,:) = exp(1i*p*x).*exp(  -((x-x0)/delta).^2  );
end
% add boundary conditions
psi(:,nx) = 0;
psi(:,1) = 0;

% set up potential
v = zeros(nx, 1);
if vtype == 1
    xmin = vpar(1);
    xmax = vpar(2);
    Vc   = vpar(3);
    v(x>xmin & x<xmax) = Vc;
    %dm(x>xmin & x<xmax) = dm(x>xmin & x<xmax) - 1/2*Vc;
end

% initialize diagonals for tridiagonal matrix A
du = zeros(nx,1); % upper diagonal  "+"
dm = zeros(nx,1); % middle diagonal "0"
dl = zeros(nx,1); % lower diagonal  "-"
f  = zeros(nx,1); % RHS

% add values to tridiagonals
du = 1/2 * 1/(dx)^2 * ones(nx,1);
dm = (1i/dt - 1/(dx)^2) * ones(nx,1) - 1/2*v;
dl = du;

%hold on;
%plot(x,v)

% fix tridiagonal boundary cases
du(2)    = 0;
dm(1)    = 1;
dm(nx)   = 1;
dl(nx-1) = 0;

% define sparse matrix:
A = spdiags([dl dm du], -1:1, nx, nx);

% iterating through n using Crank-Nichelson scheme
% some constants for simplifying f
C1 = -1/2 * 1/(dx)^2;
C2 = 1i/dt + 1/dx^2;
C3 = C1;
for n=1:nt-1
    if vtype == 1
        %disp(size(v(2:nx-1)))
        %disp(size(psi(n,2:nx-1)))
        f(2:nx-1) = C1*psi(n,3:nx) + C2*psi(n,2:nx-1) + 1/2*transpose(v(2:nx-1)).*psi(n,2:nx-1) + C3*psi(n,1:nx-2);        
    else                                                        %  vvvvv?
        f(2:nx-1) = C1*psi(n,3:nx) + C2*psi(n,2:nx-1) + C3*psi(n,1:nx-2);
    end
    psi(n+1,:) = A \ f;
    % reinforce boundary conditions
    psi(:,nx) = 0;
    psi(:,1) = 0;
end
psire = real(psi);
psiim = imag(psi);
psimod = abs(psi);

%for ti=[1:nt]
%    for xi = [1:nx]
%        prob(ti,xi) = 1/2 * sum((psimod(ti,1:xi-1).^2 + psimod(ti,2:xi).^2).*(x(2:xi) - x(1:xi-1)));
%    end
%end

prob = cumtrapz(psimod.^2, 2);
prob = prob ./ prob(:,nx);
end

