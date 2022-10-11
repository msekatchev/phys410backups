function v = totalenergy(r)
    nc = size(r,1);
    c_indecies = nchoosek(1:nc,2); % list all unique combinations of pairs of charges
        % https://www.mathworks.com/help/matlab/ref/nchoosek.html
    v = sum(vecnorm(r(c_indecies(:,1),:,:) - r(c_indecies(:,2),:,:),2,2).^-1);
end
