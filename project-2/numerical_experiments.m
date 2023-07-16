% Script for generating visualizations of 2D Schrodinger equation solutions
% for different cases.

% select a cse to generate here, from [-1,0,1,2,3,4,5,6]. Explained in
% detail in the writeup.
% Cases -1 and 0 were only used to produce example plots of the potential.
selected_case = 6;

if selected_case == -1
    idtype = 0;
    vtype  = 1;
    idpar  = [2, 3];
    vpar   = [0.5,0.7,0.2,0.5,1e6];
    tmax   = 0.05;
    lambda = 0.05;
    level = 6;

    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    
    plot_2d_potential(x,y,v,vtype,"[0.5,0.7,0.2,0.5,1e6]")
end

if selected_case == 0
    % Case 3
    idtype = 0;
    vtype  = 2;
    idpar  = [2, 3];
    vpar   = [0.2,0.4,0.6,0.9,1e6];
    tmax   = 0.05;
    lambda = 0.05;
    level = 6;

    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

    plot_2d_potential(x,y,v,vtype,"[0.2,0.4,0.6,0.9,1e6]")
end

if selected_case == 1
    % Case 1
    idtype = 0;
    vtype  = 0;
    idpar  = [2, 3];
    vpar   = 0;
    tmax   = 0.05;
    lambda = 0.05;
    level = 7;

    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    plot_2d_psire(x,y,t,psire,"Case 1")
end

if selected_case == 2
    % Case 2
    idtype = 1;
    vtype  = 0;
    idpar  = [0.5,0.5,0.075,0.075,0,0];
    vpar   = 0;
    tmax   = 0.05;
    lambda = 0.05;
    level = 8;

    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    plot_2d_psire(x,y,t,psire,"Case 2")
end

if selected_case == 3
    % Case 3
    close all;
    idtype = 1;
    vtype  = 1;
    idpar  = [0.8,0.5,0.050,0.050,-25,0];
    vpar   = [0.2,0.4,0.4,0.6,1e15];
    tmax   = 0.015;
    lambda = 0.02;
    level = 9;

    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    plot_2d_psire_potential(x,y,t,psire,"Case 3", vpar)
    %plot_2d_psire(x,y,t,psire,"Case 3")
    figure
    plot_2d_potential(x,y,v,vtype,"[0.2,0.4,0.4,0.6,1e10]")
end


if selected_case == 4
    % Case 4
    close all;
    idtype = 1;
    vtype  = 1;
    idpar  = [0.65,0.5,0.075,0.075,0,0];
    vpar   = [0.45,0.85,0.25,0.75,-1e6];
    tmax   = 0.0075;
    lambda = 0.01;
    level = 9;
    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    plot_2d_psire_potential(x,y,t,psire,"Case 4", vpar)
    %plot_2d_psire(x,y,t,psire,"Case 3")
    figure
    plot_2d_potential(x,y,v,vtype,"[0.45,0.85,0.25,0.75,-1e6]")
end

if selected_case == 5
    % Case 5
    close all;
    idtype = 1;
    vtype  = 1;
    idpar  = [0.2,0.5,0.075,0.075,20,0];
    vpar   = [0.45,0.85,0.25,0.75,-1e4];
    tmax   = 0.015;
    lambda = 0.005;
    level = 8;
    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    plot_3d_psire(x,y,t,psire,"Case 5")
    %plot_2d_psire(x,y,t,psire,"Case 3")
    figure
    plot_2d_potential(x,y,v,vtype,"[0.45,0.85,0.25,0.75,-1e4]")
end


if selected_case == 6
    % Case 6
    close all;
    idtype = 1;
    vtype  = 2;
    idpar  = [0.1,0.5,0.075,0.075,40,0];
    vpar   = [0.30,0.45,0.55,0.70,1e6];
    tmax   = 0.02;
    lambda = 0.005; 
    level = 8;
    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    plot_3d_psire(x,y,t,psire,"Case 6")
    %plot_2d_psire(x,y,t,psire,"Case 3")
    figure
    plot_2d_potential(x,y,v,vtype,"[0.30,0.45,0.55,0.70,1e4]")
end











