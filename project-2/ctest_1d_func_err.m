function [dpsi0rms, t0] = ctest_1d_func_err(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Calculate the l-2 norm of the difference of psi with the exact solution,
% to perform convergence testing of sch_1d_cn()
%
% Inputs
%
% tmax: Maximum integration time
% lmin: Minimum Discretization level
% lmax: Maximum Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parameters
%
% Outputs
%
% dpsi0rms: l-2 norm (rms) of the difference in psi between two adjacent
%           levels [nt x 1]
% t0:       Vector of time coordinates [nt]

    [x0, t0, psi0, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, level  , lambda, idtype, idpar, vtype, vpar); 

    m = idpar(1);
    
    psi_exact = sch_1d_exact(x0, t0, m);
    
    dpsi0 = psi_exact - psi0;
   
    nx = length(x0);
    
    dpsi0mod = abs(dpsi0);

    dpsi0rms = sqrt( sum(  dpsi0mod.^2  ,2) / nx ); 

end
