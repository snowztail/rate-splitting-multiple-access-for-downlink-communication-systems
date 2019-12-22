function [rate] = rs_rate(weight, bcChannel, snr, tolerance, rsRatio)
%RS_RATE Summary of this function goes here
%   Detailed explanation goes here
% MISO

[rx, tx, user] = size(bcChannel);
% reshape BC channel matrix [H] (tx * rx * user) with Ref 1
bcChannel = squeeze(permute(bcChannel, [2, 1, 3]));

% initialize common precoder using MRT & SVD (see Ref 2)
[u, ~, ~] = svd(bcChannel);
% largest left singular vector of channel matrix as common channel
comChannel = u(:, 1);
% common precoder (tx * 1)
comPrecoder = sqrt(snr * (1 - rsRatio)) * comChannel / norm(comChannel);
% private precoder (tx * user)
priPrecoder = sqrt(snr * (rsRatio / user)) * bcChannel ./ vecnorm(bcChannel);
clearvars u;

isConverged = false;
wsr = 0;
while (~isConverged)
    % compute equalizers and weights for successive precoder optimization
    [comEqualizer, priEqualizer, comWeight, priWeight, ~, ~] = rs_terms(bcChannel, comPrecoder, priPrecoder);
    % optimize common and private precoders
    [comPrecoder, priPrecoder, wsr_] = rs_solver(weight, bcChannel, snr, comEqualizer, priEqualizer, comWeight, priWeight);
    if (wsr_ - wsr) / wsr_ <= tolerance
        isConverged = true;
    end
    wsr = wsr_;
end

% compute common and private rates
[~, ~, ~, ~, comRate, priRate] = rs_terms(bcChannel, comPrecoder, priPrecoder);
% allocate common rate to different users (assume all to one user; the result is columnwise)
rate = priRate + diag(repmat(comRate, [user, 1]));

end
