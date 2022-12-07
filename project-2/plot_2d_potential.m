function plot_2d_potential(x,y,v, vtype, vpar_title)
[X,Y] = meshgrid(x,y);
contourf(X,Y,v);
xlabel("$x$",'Interpreter','latex')
ylabel("$y$",'Interpreter','latex')

title("$V(x,y)$ for vtype="+string(vtype)+", vpar = "+vpar_title,'Interpreter','latex')
colorbar
end
