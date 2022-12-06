
[x0, y0, t0, psi0, ~, ~, ~, ~] = sch_2d_adi(tmax, 2  , lambda, idtype, idpar, vtype, vpar); 

[~ ,  ~,  ~, psi1, ~, ~, ~, ~] = sch_2d_adi(tmax, 2+1, lambda, idtype, idpar, vtype, vpar);

dpsi0 = psi1(1:2:end,1:2:end,1:2:end) - psi0;
     
nx = length(x0);

dpsi0mod = abs(dpsi0);

rms1 = sqrt(mean(dpsi0mod.^2,2));

dpsi0rms = sqrt(mean(rms1.^2,3));

%dpsi0rms = sqrt( sum(  dpsi0mod.^2  ,2) / nx ); 
