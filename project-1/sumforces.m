function s = sumforces(r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sumforces: calculates the sum of the separation vectors between all 
%            charges
%
% Input arguments
%
% r: Current charge positions (nc x 3 array, where nc is the number of
% charges)
%
% Output arguments
%
% s: Array of sums of all separation vectors for each charge(nc x 3 array)
% r: Positions of charges (nc x 3 x nt array)
% v: Potential vector (length nt row vector)
% v_ec: Equivalence class counts (row vector with length determined
% by equivalence class analysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nc = size(r,1);
    s = zeros(nc,3);
    for c=1:nc
        s(c,:) = sum((r(1:end~=c,:)-r(c,:)) ./ vecnorm(r(1:end~=c,:) - r(c,:),2,2).^3);
        %                  ^^ for all charges except c
end
