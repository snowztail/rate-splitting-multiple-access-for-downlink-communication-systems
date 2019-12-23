function [rate] = rs_rate(weight, bcChannel, snr, tolerance, rsRatio)
% Function:
%   - compute the achievable user rates with rate-splitting multiple access
%
% InputArg(s):
%   - weight [u] (user * 1): user weights
%   - bcChannel [H] (rx * tx * user): broadcast channel response
%   - snr [\rho]: signal-to-noise ratio, which equals transmit power since unit noise power assumed
%   - tolerance [\epsilon]: tolerance ratio for convergence
%   - rsRatio [\alpha]: power ratio for private message streams
%
% OutputArg(s):
%   - rate: achievable user rates
%
% Comment(s):
%   - for 1-layer RS on MU-MISO systems only
%   - encode (user + 1) streams, require 1-layer SIC for all users
%   - cope with any user deployment scenario (no user grouping or ordering)
%   - maximize weighted-sum rate
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 23 Dec 19

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
