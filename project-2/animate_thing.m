function animate_thing(kmax,k_skip,x,y,y_limit,pause_time)
v = VideoWriter('psi_evolution.avi');
v.FrameRate=1;
open(v)
for k=1:k_skip:kmax
    plot(x,y(:,k))
    ylim(y_limit)
    xlim([0,1])
    frame=getframe(gcf);
    writeVideo(v,frame);
    pause(pause_time);
end
close(v)
end