function animate_thing_2df(k_skip,X,Y,Z,pause_time)
kmax = 5; % kmax = nt; each column is a timestep for me
c_limit = [min(min(min(Z))),max(max(max(Z)))];
v = VideoWriter('psi_evolution.avi');
v.FrameRate=round(1/pause_time);
open(v)
for k=1:k_skip:kmax
    %size_x = size(X)
    %size_y = size(Y)
    %size_z = size(permute(Z(k,:,:),[3 2 1]))
    pcolor(X,Y,squeeze(Z(k,:,:)))
    shading interp;
    colorbar
    %grid
    %size(x)
    %size(y(:,k))
    %11
    xlim([0,1])
    xlim([0,1])
    %colorbar.clim(c_limit)
    frame=getframe(gcf);
    writeVideo(v,frame);
    pause(pause_time);
end
close(v)
end