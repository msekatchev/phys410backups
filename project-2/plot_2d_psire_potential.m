function plot_2d_psire_potential(x,y,t,psire,case_name,vpar)
% Create a video of the evolution of the solution on a 2d contourf plot

x1 = vpar(1);
x2 = vpar(2);
y1 = vpar(3);
y2 = vpar(4);

box_y = y1:0.01:y2;
box_x = x1:0.01:x2;
len_y = length(box_y);
len_x = length(box_x);
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
        plot(x1*ones(len_y,1),box_y,'k-')
        plot(x2*ones(len_y,1),box_y,'k-')
        plot(box_x,y1*ones(len_x,1),'k-')
        plot(box_x,y2*ones(len_x,1),'k-')
        title(sprintf('Time step %d of %d: t = %.3g', ti, nt, t(ti)));
        for frame = 1:50
            writeVideo(vid, getframe(gcf));
        end
    end


    
    contourf(X,Y,squeeze(psire(ti,:,:)));
    title(sprintf('Time step %d of %d: t = %.3g', ti, nt, t(ti)));
    
    plot(x1*ones(len_y,1),box_y,'k-')
    plot(x2*ones(len_y,1),box_y,'k-')
    plot(box_x,y1*ones(len_x,1),'k-')
    plot(box_x,y2*ones(len_x,1),'k-')
    drawnow;
    writeVideo(vid,getframe(gcf));
end
close(vid);
end
