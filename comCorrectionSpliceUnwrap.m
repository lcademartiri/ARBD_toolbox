function [spt,sput]=comCorrectionSpliceUnwrap(p,bins,S)
	fprintf('### Initializing: COM-Correction and Splicing ###\n');
    T_steps=size(p,3);
    COMw = mean(p,1);   % wrapped COM
    pt=p-COMw; % COM-subtracted wrapped coordinates
    slices=(1:T_steps/bins:T_steps)';
    for i0=1:slices
        sp{i0,1}=p(:,:,slices(i0):slices(i0)+T_steps/bins-1); % sliced non-COM-subtracted wrapped trajectories
        spt{i0,1}=pt(:,:,slices(i0):slices(i0)+T_steps/bins-1); % sliced COM-subtracted wrapped trajectories
    end
    if S.bc==2
        for irep=1:bins
            spu{irep,1} = zeros(size(sp{irep,1})); % initialize unwrapped positions array
            spu{irep,1}(:,:,1) = sp{irep,1}(:,:,1); % starting positions
            for t = 2:T_steps/10
                dp = sp{irep,1}(:,:,t) - sp{irep,1}(:,:,t-1); % displacements
                dp(dp >  S.br) = dp(dp >  S.br) - 2*S.br; % MIC #1
                dp(dp < -S.br) = dp(dp < -S.br) + 2*S.br; % MIC #1
                spu{irep,1}(:,:,t) = spu{irep,1}(:,:,t-1) + dp; % sliced non-COM-subtracted UNwrapped trajectories
            end
            COMu = (mean(spu{irep,1},1));
            sput{irep,1}=spu{irep,1}-COMu; % sliced COM-subtracted UNwrapped trajectories
        end
    elseif S.bc==3
        [~,rotmat]=FCCrotate([1,0,0],[1,1,1]./norm([1,1,1]));
        rotmat=rotmat';
        for irep=1:bins
            spu{irep,1} = zeros(size(sp{irep,1})); % initialize unwrapped positions array
            spu{irep,1}(:,:,1) = sp{irep,1}(:,:,1); % starting positions
            A_mat = S.fcc.A; 
            invA_mat = S.fcc.invA;
            for t = 2:T_steps/10
                dp = sp{irep,1}(:,:,t) - sp{irep,1}(:,:,t-1); % A. Calculate raw Cartesian displacement
                dp_frac = dp * invA_mat; % B. Convert displacement to Fractional coordinates
                dp_frac = dp_frac - round(dp_frac); % C. Apply Minimum Image Convention in Fractional Space
                dp_mic = dp_frac * A_mat; % D. Convert back to Cartesian displacement
                spu{irep,1}(:,:,t) = spu{irep,1}(:,:,t-1) + dp_mic; % E. Accumulate
            end
    
            COMu = mean(spu{irep,1}, 1);
            sput{irep,1} = spu{irep,1} - COMu;
        end
    end
	fprintf('=== Completed: COM-Correction and Splicing === \n');
end