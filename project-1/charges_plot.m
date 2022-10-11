function charges_plot(t, r, interactive)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charges_plot(t, r, interactive) plots charge positions.  If interactive
% is non-zero, waits for user input after each time step.
% 
% Input arguments
%
%  t:    Vector of output times (length nt)
%  r:    Charge positions (3d array of size nc x 3 x nt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [nc, ~, nt] = size(r);
   for it = 1 : nt
      figure(1);
      clf;
      x = squeeze(r(:,1,it));
      y = squeeze(r(:,2,it));
      z = squeeze(r(:,3,it));
      plot3(x, y, z, '.r', 'MarkerSize', 20);
      pbaspect([1,1,1]);
      axis([-1 1 -1 1 -1 1]);
      title(sprintf('Time step %d of %d: t = %.3g', it, nt, t(it)));
      drawnow;
      if interactive
         dontcare = input('Enter anything to continue: ', 's');
      end
   end
end
