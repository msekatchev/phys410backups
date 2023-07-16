function plot_2d_psire(x,y,t,psire,case_name)
% Create a video of the evolution of the solution on a 2d contourf plot
%
% Inputs:
% x: Vector of x coordinates [nx]
% y: Vector of y coordinates [ny]
% t: Vector of t coordinates [nt]
% psire Array of computed psi_re values [nt x nx x ny]
% case_name: string for the name of the video
% 
% Outputs:
% Saves a video titled "case_name.avi", with the evolution of psire.

figure
hold on
xlabel("$x$",'Interpreter','latex')
ylabel("$y$",'Interpreter','latex')

nt = length(t);
max_psire = max(psire,[],'all');
min_psire = min(psire,[],'all');
% Initialize and open the video object ...
vid = VideoWriter(case_name+'.avi');
open(vid);

[X,Y] = meshgrid(x,y);

for ti = 1 : 1 : nt
    figure(1);
    
    disp(ti)
    if ti==1
        xlabel("$x$",'Interpreter','latex')
        ylabel("$y$",'Interpreter','latex')
        colorbar
        colormap("turbo")
        caxis([min_psire,max_psire])  
        % display initial solution for a few seconds
        contourf(X,Y,squeeze(psire(ti,:,:)));
        title(sprintf('Time step %d of %d: t = %.3g', ti, nt, t(ti)));
        for frame = 1:50
            writeVideo(vid, getframe(gcf));
        end
    end


    
    contourf(X,Y,squeeze(psire(ti,:,:)));
    title(sprintf('Time step %d of %d: t = %.3g', ti, nt, t(ti)));

    
    drawnow;
    writeVideo(vid,getframe(gcf));
end
close(vid);
end
