function r0 = random_r0(nc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random_r0: Generates random charge distribution constrained on sphere
%
% Input arguments
%
% nc: Number of charges
%
% Output arguments
%
% r0: Array of initial random positions of nc charges (nc x 3 array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r0 = 2*rand(nc,3)-1;
    r0 = r0 ./ sqrt(r0(:,1).^2 + r0(:,2).^2 + r0(:,3).^2);
    % implementation of the above line tested using
    %   plot(r0(:,1),r0(:,2),"o") before and after its addition.
end
    