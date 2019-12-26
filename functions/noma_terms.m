function [equalizer, mmseWeight, rate] = noma_terms(bcChannel, precoder, order)
% Function:
%   - construct power terms, equalizers, MMSEs, and weights for precoder optimization
%   - calculate achievable user rates
%
% InputArg(s):
%   - bcChannel [H] (rx * tx * user): broadcast channel response
%   - precoder [p] (tx * user): optimum precoders maximizing WSR
%   - order [\pi] (user * 1): decoding sequence
%
% OutputArg(s):
%   - equalizer [g] (user * 1): optimum MMSE equalizer
%   - mmseWeight [u^mmse] (user * 1): optimum MMSE weights
%   - rate [R_i] (user * 1): achievable user rates corresponding to maximum WSR
%
% Comment(s):
%   - SIC stops after users decoding own layer
%   - the rate of each layer must be achievable for those decode it (thus the minimum of all user rates on this layer)
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 25 Dec 19


[~, user] = size(bcChannel);
% receive power terms [T]
powTerm = zeros(user);
% interference power terms [I]
intPowTerm = zeros(user);

% total receive power [T^c] (i.e. power at the first layer)
powTerm(:, 1) = sum(abs(bcChannel' * precoder) .^ 2, 2) + 1;
for iUser = 1 : user
    for iLayer = 1 : user
        if iLayer == 1
            % interference power at the first layer
            intPowTerm(iUser, iLayer) = powTerm(iUser, iLayer) - abs(bcChannel(:, iUser)' * precoder(:, order(iLayer))) .^ 2;
        else
            % remaining power at the i-th layer
            powTerm(iUser, iLayer) = powTerm(iUser, iLayer - 1) - abs(bcChannel(:, iUser)' * precoder(:, order(iLayer - 1))) .^ 2;
            % interference power at the i-th layer
            intPowTerm(iUser, iLayer) = intPowTerm(iUser, iLayer - 1) - abs(bcChannel(:, iUser)' * precoder(:, order(iLayer))) .^ 2;
        end
    end
end

% SIC stops after users decoding own layer
powTerm(order, :) = tril(powTerm(order, :));
powTerm(powTerm == 0) = NaN;
intPowTerm(order, :) = tril(intPowTerm(order, :));
intPowTerm(intPowTerm == 0) = NaN;

% optimum MMSE equalizers [g]
equalizer = zeros(user);
for iUser = 1 : user
    for iLayer = 1 : user
        equalizer(iUser, iLayer) = precoder(:, iLayer)' * bcChannel(:, iUser) / powTerm(iUser, iLayer);
    end
end

% MMSEs [\epsilon]
mmse = powTerm .\ intPowTerm;

% optimum MMSE weights
mmseWeight = 1 ./ mmse;

% corresponding maximum achievable rates at each layer
rate = real(log2(powTerm ./ intPowTerm));
% the rate of each layer must be achievable for those decode it (thus the minimum of all user rates on this layer)
rate = min(rate);

end
