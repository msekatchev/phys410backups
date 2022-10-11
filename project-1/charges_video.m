function charges_video(t, r, filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charges_video(t, r) plots charge positions and makes an 
% AVI video file named 'charges.avi'.
% 
% Input arguments
%
%  t:        Vector of output times (length nt)
%  r:        Charge positions (3d array of size nc x 3 x nt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get the number of time steps from the size of the coordinate array ...
   [~, ~, nt] = size(r);

   % Initialize and open the video object ...
   v = VideoWriter('charges.avi');
   open(v);

   % For each time step ...
   for it = 1 : nt
      % Make plot of charges ...
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
      if it == 1 
         % Record 5 seconds of the first frame ... 
         for frame = 1 : 75
            writeVideo(v, getframe(gcf));
         end
      else
         % Record the frame ...
         writeVideo(v, getframe(gcf));
      end
   end
   % Close the video object ...
   close(v);
end
