function [mask,N_eff]=sphericalMask(S,spt)
    % Mask is defined on the COM-corrected 'p'
    if S.bc == 2 % PBC
        for irep=1:size(spt,1)
            dist_sq = sum(spt{irep,1}.^2, 2); 
            mask{irep,1} = squeeze(dist_sq < S.br^2); 
            N_eff{irep,1} = mean(sum(mask{irep,1}, 1));
        end
    elseif S.bc==3
        for irep=1:size(spt,1)
            box_heights = abs(dot(S.fcc.A, S.fcc.normals, 2));
            R_mask = min(box_heights) / 2; % The inscribed radius is half the smallest height
            dist_sq = sum(spt{irep,1}.^2, 2); % (N x 1 x T_seg)
            mask{irep,1} = squeeze(dist_sq < R_mask^2);
            N_eff_t = sum(mask{irep,1}, 1);
            N_eff{irep,1} = mean(N_eff_t);
            if N_eff{irep,1} < 1, N_eff{irep,1} = 1; end
        end
    elseif S.bc==1 % SBC
        for irep=1:size(spt,1)
            mask{irep,1} = true(S.N, 1, T_steps/10);
            N_eff{irep,1} = S.N;
        end
    end
end