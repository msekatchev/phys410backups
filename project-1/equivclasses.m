function v_ec = equivclasses(r, epsec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equivclasses: produce a vector defining the number of charges in each
%               equivalence class.
%
% Input arguments
%
% r: Final charge positions (nc x 3 array, where nc is the number of
% charges)
% epsec: Tolerance for this equivalence class analysis
%
% Output arguments
%
% v_ec: vector defining number of charges in each equivalence class
%       (1 x equiv array, where equiv is the number of equivalence classes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nc = size(r,1);
    v_ec = zeros(1,nc);
    d = zeros(nc,nc,3); % create nc x nc x 3 matrix for the difference calculations
    for c=1:nc
        d(c,:,:) = ones(nc,1)*r(c,:); % populate matrix with charge positions. The 
                                      % 3rd dimension allows to store x, y,
                                      % z values of charges separately.
    end
    d = vecnorm(d-pagetranspose(d),2,3); % calculate distances detween all charges.
                                         % this is d_ij = |r_j - r_i| from
                                         % writeup
    
    dsorted = sort(d,2);                 % sort distances array in ascending order
    matched_charges = ones(nc,1); % simple list of charges, m(i) will be replaced with 
                                  % 0 if the charge has already been used
                                  % in an equivalence class
    for c=1:nc
        % for each charge, calculate the abs of the difference in sorted distances
        % between all other charges and check for equivalences
        equivalence_matrix = abs(dsorted - ones(nc,nc).*dsorted(c,:));
        
        for m=1:nc
            % check if charge has not been used in equivalence class
            if matched_charges(m) ~= 0 && all(equivalence_matrix(m,:) <= epsec)
            % and check if |dbar_i-dbar_i'| <= epsec element-wise for all
            % sorted differences
                v_ec(c) = v_ec(c) + 1;
                matched_charges(m) = 0;
            end
        end
    end
    v_ec = sort(v_ec,'descend');
end