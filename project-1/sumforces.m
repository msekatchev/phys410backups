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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nc = size(r,1);
    s = zeros(nc,3);
    for c=1:nc
        s(c,:) = sum((r(1:end~=c,:)-r(c,:)) ./ vecnorm(r(1:end~=c,:) - r(c,:),2,2).^3);
        %                  ^^ for all charges except c
end
