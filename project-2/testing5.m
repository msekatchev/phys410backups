tmax = 0.05;
level = 9;
lambda = 0.05;
idtype = 0;
idpar = [2,3];

vpar = 0;
vtype = 0;

% initialize discretization mesh
nx = 2^level + 1;
ny = nx;
x = linspace(0,1,nx);
y = linspace(0,1,ny);
dx = x(2) - x(1);
dy = dx;
dt = lambda * dx;
nt = round(tmax / dt) + 1;

t   = [0: nt-1] * dt;

psi = zeros(nt,nx,ny, 'like', 1j);

y_T = y.';


% initialize psi depending on IC type -------------------------------------
if idtype == 0
    mx     = idpar(1);
    my     = idpar(2);

    psi(1,:,:) = sin(mx*pi*x).*sin(my*pi*y_T);
else
    x0     = idpar(1);
    y0     = idpar(2);
    deltax = idpar(3);
    deltay = idpar(4);
    px     = idpar(5);
    py     = idpar(6);
    delta = idpar(2);
    
    psi(1,:,:) = exp(1i*px*x).*exp(1i*py*y_T).*exp(-( ((x-x0)/deltax).^2 + ((y_T-y0)/deltay).^2));
end

% add boundary conditions -------------------------------------------------
psi(:,nx, :) = 0;
psi(:, 1, :) = 0;
psi(:, :,ny) = 0;
psi(:, :, 1) = 0;

% set up potential depending on vtype -------------------------------------
v = zeros(nx, ny);
if vtype == 1
    xmin = vpar(1);
    xmax = vpar(2);
    ymin = vpar(3);
    ymax = vpar(4);
    Vc   = vpar(5);
    
    v((x>xmin & x<xmax) & (y_T>ymin & y_T<ymax)) = Vc;
end
if vtype == 2
    x1 = vpar(1);
    x2 = vpar(2);
    x3 = vpar(3);
    x4 = vpar(4);
    Vc = vpar(5);   
    
    % j' and j'+1
    j_slit0 = (ny-1)/4 + 1;
    j_slit1 = j_slit0+1;
    
    v(:,j_slit0:j_slit1) = Vc;
    v(((x1<=x) & x<=x2) | ((x3<=x) & (x<=x4)),:) = 0;
end


% LHS - A matrix ----------------------------------------------------------

% initialize diagonals for tridiagonal matrix A
du = zeros(nx,1); % upper diagonal  "+"
dm = zeros(nx,1); % middle diagonal "0"
dl = zeros(nx,1); % lower diagonal  "-"

% add values to tridiagonals
du =   - 1i * dt/(2*dx^2) * ones(nx,1);
dm = 1 + 1i * dt/(  dx^2) * ones(nx,1);
dl = du;

% fix tridiagonal boundary cases
du(2)    = 0;
dm(1)    = 1;
dm(nx)   = 1;
dl(nx-1) = 0;

% define sparse matrix:
A = spdiags([dl dm du], -1:1, nx, nx);


% RHS - G constants -------------------------------------------------------
Gp =     1i * dt/(2*dy^2) * ones(nx,ny);
G0 = 1 - 1i * dt/(  dy^2) * ones(nx,ny) - 1i*dt/2*v;
Gm = Gp;

% RHS - F constants -------------------------------------------------------
Fp =     1i * dt/(2*dx^2) * ones(nx,ny);
F0 = 1 - 1i * dt/(  dx^2) * ones(nx,ny);
Fm = Fp;

% LHS - B constants -------------------------------------------------------
Bp =   - 1i * dt/(2*dy^2) * ones(1,ny);
B0 = 1 + 1i * dt/(  dy^2) + 1i*dt/2*v;     % this one is [nx,ny], a different row will be selected for each ii
Bm = Bp;

G = zeros(nx,ny);
F = zeros(nx,ny);
psi_temp = zeros(nx,ny, 'like', 1j);

psi_nplus15 = zeros(nt,nx,ny, 'like', 1j);

for n = [1:nt-1]
    psi_n = squeeze(psi(n,:,:));
    G(:,(2:ny-1)) = Gp(:,(3:ny)) .* psi_n(:,(3:ny)) + G0(:,(2:ny-1)) .* psi_n(:,(2:ny-1)) + Gm(:,(1:ny-2)) .* psi_n(:,(1:ny-2));

    F((2:nx-1),:) = Fp((3:nx),:) .* G((3:nx),:) + F0((2:nx-1),:) .* G((2:nx-1),:) + Fm((1:nx-2),:) .* G((1:nx-2),:);
    
    
    for jj = [2:ny-1]
        psi_temp(:,jj) = A \ F(:,jj);
    end
    
    % reinforce BCs
    psi_temp(nx, :) = 0;
    psi_temp( 1, :) = 0;
    psi_temp( :,ny) = 0;
    psi_temp( :, 1) = 0;
    
    for ii = [2:nx-1]
        % define sparse matrix for B, different one for each ii
        Bupper = Bp.';
        Bmain  = B0(ii,:).';
        Blower = Bm.';
        % fix tridiagonal boundary cases
        Bupper(2)    = 0;
        Bmain(1)     = 1;
        Bmain(nx)    = 1;
        Blower(nx-1) = 0; 
        
        B = spdiags([Blower Bmain Bupper], -1:1, ny, ny);
        
        psi(n+1,ii,:) = (B \ psi_temp(ii,:).').';
    end
    
    % reinforce BCs
    psi(:,nx, :) = 0;
    psi(:, 1, :) = 0;
    psi(:, :,ny) = 0;
    psi(:, :, 1) = 0;
    
end


probability = zeros(1,nt);

for n = 1:nt
    probability(n) = norm(squeeze(psi(n,:,:)));
end

psire = real(psi);
psiim = imag(psi);
psimod = abs(psi);

figure;
hold on
plot(x,psire(20,:,20));
psi_exact = sch_2d_exact(x,y,t,mx,my);
plot(x,psi_exact(20,:,20));
hold off;
figure;
plot(t,probability);
