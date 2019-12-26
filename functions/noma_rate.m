function [rate] = noma_rate(weight, bcChannel, snr, tolerance)
% Function:
%   - compute the achievable user rates for non-orthogonal multiple access with successive interference cancellation
%
% InputArg(s):
%   - weight [u] (user * 1): user weights
%   - bcChannel [H] (rx * tx * user): broadcast channel response
%   - snr [\rho]: signal-to-noise ratio, which equals transmit power since unit noise power assumed
%   - tolerance [\epsilon]: tolerance ratio for convergence
%
% OutputArg(s):
%   - rate: achievable user rates
%
% Comment(s):
%   - as a special case of RSMA
%   - for NOMA with ordered-SIC on MU-MISO systems only
%   - encode (user) streams, subject to error propagation
%   - require aligned users
%   - examine all possible user ordering for optimal performance
%   - maximize weighted-sum rate
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 25 Dec 19


[rx, tx, user] = size(bcChannel);
% reshape BC channel matrix [H] (tx * rx * user) with Ref 1
bcChannel = squeeze(permute(bcChannel, [2, 1, 3]));

% initialzie precoders (tx * user)
precoder = sqrt(snr / user) * bcChannel ./ vecnorm(bcChannel);

% all possible user orders [\pi] (permutations * user)
order = perms(1 : user);
% number of permutations
nPerms = size(order, 1);

rate = zeros(nPerms, user);
for iPerm = 1 : nPerms
    isConverged = false;
    wsr = 0;
    while(~isConverged)
        % compute equalizers and weights for successive precoder optimization
        [equalizer, mmseWeight, ~] = noma_terms(bcChannel, precoder, order(iPerm, :));
        % optimize precoders
        [precoder, wsr_] = noma_solver(weight, bcChannel, snr, equalizer, mmseWeight, order(iPerm, :));
        if (wsr_ - wsr) / wsr_ <= tolerance
            isConverged = true;
        end
        wsr = wsr_;
    end
    % compute achievable user rates for the corresponding decoding order
    [~, ~, rate(iPerm, :)] = noma_terms(bcChannel, precoder, order(iPerm, :));
end

end
