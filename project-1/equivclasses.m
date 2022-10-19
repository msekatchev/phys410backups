function v_ec = equivclasses(r, epsec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equivclasses: XXXX
%
% Input arguments
%
% r: Final charge positions (nc x 3 array, where nc is the number of
% charges)
%
% Output arguments
%
% v_ec: XXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nc = size(r,1);
    v_ec = zeros(1,nc);
    d = zeros(nc,nc,3);
    for c=1:nc
        d(c,:,:) = ones(nc,1)*r(c,:);
    end
    d = vecnorm(d-pagetranspose(d),2,3);
    dsorted = sort(d,2);
    matched_charges = ones(nc,1);

    for c=1:nc
        difference_vector = vecnorm(dsorted - ones(nc,nc).*dsorted(c,:),2,2)
        for m=1:nc
            if matched_charges(m) ~= 0
                if difference_vector(m) <= epsec
                    v_ec(c) = v_ec(c) + 1
                    matched_charges(m) = 0
                end
            end
        end
            
    end
    
%    pairs = nchoosek(1:nc,2);
%    for pair1 = 1:length(pairs)
%        for pair2 = 1:length(pairs)
%            dsorted(pairs(pair1,))
%    end
end