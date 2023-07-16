function plot_2d_potential(x,y,v, vtype, vpar_title)
% Create a 2D contourf plot of the potential
%
% Inputs:
% x: Vector of x coordinates [nx]
% y: Vector of y coordinates [ny]
% v: Array of potential values [nx x ny]
% psire: Array of computed psi_re values [nt x nx x ny]
% vpar_title: string with the vpar parameters, used in the title of the
%   plot
% 
[X,Y] = meshgrid(x,y);
contourf(X,Y,v);
xlabel("$x$",'Interpreter','latex')
ylabel("$y$",'Interpreter','latex')

title("$V(x,y)$ for vtype="+string(vtype)+", vpar = "+vpar_title,'Interpreter','latex')
colorbar
end
