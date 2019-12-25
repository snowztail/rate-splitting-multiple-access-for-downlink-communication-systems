function [rate] = noma_rate(weight, bcChannel, snr, tolerance)
%NOMA_RATE Summary of this function goes here
%   Detailed explanation goes here


[rx, tx, user] = size(bcChannel);
% reshape BC channel matrix [H] (tx * rx * user) with Ref 1
bcChannel = squeeze(permute(bcChannel, [2, 1, 3]));

% initialzie precoders (tx * user)
precoder = sqrt(snr / user) * bcChannel ./ vecnorm(bcChannel);

% all possible user orders [\pi] (permutations * user)
order = perms(1 : user);
% number of permutations
nPerms = size(order, 1);

for iPerm = 1 : nPerms
    isConverged = false;
    wsr = 0;
    while(~isConverged)
    % compute equalizers and weights for successive precoder optimization
    [equalizer, mmseWeight, ~] = noma_terms(bcChannel, precoder, [1 2]);
    % optimize precoders
    [precoder, wsr] = noma_solver(weight, bcChannel, snr, equalizer, mmseWeight, [1 2])

    end
end









% isConverged = false;
% wsr = 0;
% while (~isConverged)
%     % compute equalizers and weights for successive precoder optimization
%     [comEqualizer, priEqualizer, comWeight, priWeight, ~, ~] = rs_terms(bcChannel, comPrecoder, precoder);
%     % optimize common and private precoders
%     [comPrecoder, precoder, wsr_] = rs_solver(weight, bcChannel, snr, comEqualizer, priEqualizer, comWeight, priWeight);
%     if (wsr_ - wsr) / wsr_ <= tolerance
%         isConverged = true;
%     end
%     wsr = wsr_;
% end

% % compute common and private rates
% [~, ~, ~, ~, comRate, priRate] = rs_terms(bcChannel, comPrecoder, precoder);
% % allocate common rate to different users (assume all to one user; the result is columnwise)
% rate = priRate + diag(repmat(comRate, [user, 1]));
end
