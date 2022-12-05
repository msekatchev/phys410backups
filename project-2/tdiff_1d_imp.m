% tdiff_1d_imp.m: Driver for diff_1d_imp ... Solution of 1d diffusion
% equation using O(dt, dx^2) implicit scheme.
more off;

format long;

% idtype = 0   ->  sine initial data (exact solution known)
% idtype = 1   ->  gaussian initial data 

idtype = 0 

tmax = 0.2
x0 = 0.5
delta = 0.05
omega = 2 * pi

minlevel = 6
maxlevel = 9

olevel = 6
ofreq  = 1

dtbydx = 0.1

% Enable for plotting ...
plotit = 1
if plotit
   close all
end

% Perform computation at various levels of discretization, store
% results in cell arrays ...
for l = minlevel : maxlevel
   tstart = tic;
   % Compute the solution  ...
   [x{l}, t{l}, u{l}] = diff_1d_imp(tmax, l, dtbydx, omega, ...
                                    x0, delta, idtype, 0);
   telapsed = toc(tstart)
   [nt{l}, nx{l}] = size(u{l});
   stride{l} = ofreq * 2^(l - olevel);
   
   % If possible, compute exact solution ...
   if idtype == 0
      uxct{l} = zeros(nt{l}, nx{l});
      for n = 1 : nt{l}
         uxct{l}(n,:) = exp(-omega^2 * t{l}(n)) * sin(omega * x{l});
      end
   end

   if plotit
      for it = 1 : nt{l}
         figure(l);
         plot(x{l},u{l}(it,:));
         xlim([0, 1]);
         ylim([-1, 1]);
         drawnow;
      end
   end
   % If possible, compute and output error norm ...
   if idtype == 0
      uerr = u{l} - uxct{l};
      fprintf('Level: %d  ||error||: %25.16e\n', ...
         l, sqrt(sum(sum((uerr .* uerr))) / prod(size(uerr))));
   end
end
