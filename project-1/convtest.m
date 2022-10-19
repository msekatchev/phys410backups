nc = 4;
tmax = 10;
gamma = 1;
epsec = 1.0e-5;

r0 = [ [1, 0, 0]; [0, 1, 0]; [0, 0, 1]; (sqrt(3)/3) * [1, 1, 1] ];

levels = [10, 11, 12, 13];

[t10, r10, ~, ~] = charges(r0, tmax, 10, gamma, epsec);
[t11, r11, ~, ~] = charges(r0, tmax, 11, gamma, epsec);
[t12, r12, ~, ~] = charges(r0, tmax, 12, gamma, epsec);
[t13, r13, ~, ~] = charges(r0, tmax, 13, gamma, epsec);
x10 = r10(1,1,:); x11 = r11(1,1,:); x12 = r12(1,1,:); x13 = r13(1,1,:);











%for l = 10:13
%    nt = 2^l + 1;
%    [t, r, ~, ~] = charges(r0, tmax, level, gamma, epsec);
%    x_fine = r(1,1,:);
%    [t, r, ~, ~] = charges(r0, tmax, level, gamma, epsec);
%    x_coarse = 
%    %x(1:2:end)
%    
%    
%    hold on;
%    plot(x, t);
%    
%    x = x(1:2:end)
%    t = t(1:2:end)
%    
%end