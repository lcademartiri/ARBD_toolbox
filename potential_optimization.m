clearvars -except S H*
[~, cmdout] = system('wmic cpu get L2CacheSize, L3CacheSize /value');
tokens = regexp(cmdout, '\d+', 'match');
cacheSizeMB = max(str2double(tokens))/1024;
for mn=1:14
    
    p=rand(2^mn,3)*S.br*2-S.br;
    N(mn,1)=2^mn;
    for irep=1:5
        tic
            [disppot2,pairs_i,pairs_j,d_mic_out] = potential_displacements_v2(p, S, H,H_interpolant, 0);
        times2(mn,irep)=toc;
    end
    for irep=1:5    
        tic
            [disppot3,pairs_i,pairs_j,d_mic_out] = potential_displacements_v3(p, S, H,H_interpolant, 0);
        times3(mn,irep)=toc;
    end
    for irep=1:5    
        % tic
        %     [disppot4,pairs_i,pairs_j,d_mic_out] = potential_displacements_v4(p, S, H,H_interpolant, 0);
        % times4(mn,irep)=toc;
    end
    for irep=1:5
        % tic
        %     [disppot5,pairs_i,pairs_j,d_mic_out] = potential_displacements_v5(p, S, H,H_interpolant);
        % times5(mn,irep)=toc;
    end
    for irep=1:5
        % tic
        %     [disppot6,pairs_i,pairs_j,d_mic_out] = potential_displacements_v6(p, S, H,H_interpolant);
        % times6(mn,irep)=toc;
    end
    for irep=1:5
        tic
            [disppot7,pairs_i,pairs_j,d_mic_out] = potential_displacements_v7(p, S, H, H_interpolant, 0, cacheSizeMB);
        times7(mn,irep)=toc;
    end
        disp([mn,irep])
    
end
times(:,1)=N;
times(:,2)=mean(times2,2);
times(:,3)=mean(times3,2);
% times(:,4)=mean(times4,2);
% times(:,5)=mean(times5,2);
% times(:,6)=mean(times6,2);
times(:,7)=mean(times7,2);
figure
loglog(times(:,1),times(:,2))
hold
loglog(times(:,1),times(:,3))
% loglog(times(:,1),times(:,4))
% loglog(times(:,1),times(:,5))
% loglog(times(:,1),times(:,6))
loglog(times(:,1),times(:,7))
lgd = legend('v2', 'v3', 'v4', 'v5', 'v6', 'v7');