function Fe = excess_frac_prob(tmax, level, lambda, idtype, idpar, vtype, vpar, x1, x2)
% Compute the excess fractional probability for sch_1d_cn() between the
% specified range and with the specified parameters.
% Inputs
%
% tmax: Maximum integration time
% level: Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parameters
% x1: lower limit for the excess fractional probability calculation
% x2: uppwer limit for the excess fractional probability calculation
%
% Outputs
%
% Fe: excess fractional probability between [x1,x2]

[x, ~, ~, ~, ~, ~, prob, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

nx = length(x);
% calculate temporal average
Pbar = mean(prob);
% normalize temporal average
Pbar = Pbar / Pbar(nx);

loc_nearest_x1 = abs(x-x1)==min(abs(x-x1));
loc_nearest_x2 = abs(x-x2)==min(abs(x-x2));

nearest_x1 = x(loc_nearest_x1);
nearest_x2 = x(loc_nearest_x2);

Fe = (Pbar(loc_nearest_x2) - Pbar(loc_nearest_x1)) ./ (nearest_x2 - nearest_x1);    
end
