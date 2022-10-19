function [t, r, v, v_ec] = charges(r0, tmax, level, gamma, epsec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charges: Top-level function for solution of charges-on-a-sphere
% problem.
%
% Input arguments
%
% r0: Initial positions (nc x 3 array, where nc is the number of
% charges)
% tmax: Maximum simulation time
% level: Discretization level
% gamma: Dissipation coefficient
% epsec: Tolerance for equivalence class analysis
%
% Output arguments
%
% t: Vector of simulation times (length nt row vector)
% r: Positions of charges (nc x 3 x nt array)
% v: Potential vector (length nt row vector)
% v_ec: Equivalence class counts (row vector with length determined
% by equivalence class analysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nc = size(r0,1);
    
    nt = 2^level+1;
    dt = tmax ./ (nt - 1);
    t = 0:dt:tmax;
    
    r = zeros(nc,3,nt);
    v = zeros(length(t),1);
    
    % initial conditions r^1 = r^2
    r(:,:,1) = r0;
    r(:,:,2) = r0;
    
    % alpha and beta for simplifying & minimizing calculations
    a = gamma / (2*dt);
    b = 1/dt^2;
    % ^ these don't need to be computed anew for each iteration
    
    v(1) = totalenergy(r(:,:,1));
    
    for n=2:1:nt-1
        r(:,:,n+1) = (1/(a+b))*( -sumforces(r(:,:,n)) + 2*b*r(:,:,n) + (a-b)*r(:,:,n-1));
        r(:,:,n+1) = r(:,:,n+1) ./ vecnorm(r(:,:,n+1),2,2);
        v(n) = totalenergy(r(:,:,n));
    end
    
    v(end) = totalenergy(r(:,:,end));
    
    v_ec = equivclasses(r(:,:,end), epsec);
    %v_ec = 0;
end
