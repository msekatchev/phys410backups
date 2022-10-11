function s = sumforces(r)
    nc = size(r,1);
    s = zeros(nc,3);
    for c=1:nc
        s(c,:) = sum((r(1:end~=c,:)-r(c,:)) ./ vecnorm(r(1:end~=c,:) - r(c,:),2,2).^3);
        %                ^^ for all charges except c
end
