selected_case = 3;

if selected_case == -1
    % Case 2
    idtype = 0;
    vtype  = 1;
    idpar  = [2, 3];
    vpar   = [0.5,0.7,0.2,0.5,1e6];
    tmax   = 0.05;
    lambda = 0.05;
    level = 6;

    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    
    %plot_2d_psire(x,y,t,psire,"Case 0")
    %figure
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
    % Case 1
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
    % Case 1
    close all;
    idtype = 1;
    vtype  = 1;
    idpar  = [0.8,0.5,0.050,0.050,-25,0];
    vpar   = [0.2,0.4,0.4,0.6,1e15];
    tmax   = 0.015;
    lambda = 0.02;
    level = 9;

    [x y t psi psire psiim psimod v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    plot_2d_psire_potential(x,y,t,psire,"Case 3--", vpar)
    %plot_2d_psire(x,y,t,psire,"Case 3")
    %figure
    %plot_2d_potential(x,y,v,vtype,"[0.2,0.4,0.4,0.6,1e10]")
end













