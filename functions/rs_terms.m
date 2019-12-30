function [comEqualizer, priEqualizer, comWeight, priWeight, comRate, priRate] = rs_terms(bcChannel, comPrecoder, priPrecoder)
% Function:
%   - construct power terms, equalizers, MMSEs, and weights for precoder optimization
%   - calculate achievable common and private rates
%
% InputArg(s):
%   - bcChannel [H] (rx * tx * user): broadcast channel response
%   - comPrecoder [p_i^c] (tx * 1): optimum common precoder
%   - priPrecoder [p_i^i] (tx * user): optimum private precoders
%
% OutputArg(s):
%   - comEqualizer [g_i^c] (user * 1): optimum MMSE equalizers for common stream
%   - priEqualizer [g_i^i] (user * 1): optimum MMSE equalizers for (self) private stream
%   - comWeight [u_i^c] (user * 1): optimum MMSE weights for common stream
%   - priWeight [u_i^i] (user * 1): optimum MMSE weights for (self) private stream
%   - comRate [R_c]: achievable common rate
%   - priRate [R_i] (1 * user): user private rates
%
% Comment(s):
%   - require sorted channels and precoders based on decoding order (hence output is also sorted)
%   - user-i decodes layer-1 to layer-i (hence elements on and below diagonal are vaild)
%   - SIC stops after users decoding own layer
%   - the private rate of each layer must be achievable for those decode it
%   - the common rate must be achievable for all users
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 27 Dec 19


[~, user] = size(bcChannel);

% private stream power [T_i^i, I_i^c] ($user$ layers)
priPowTerm = zeros(user);
% total private power at the first layer
priPowTerm(:, 1) = sum(abs(bcChannel' * priPrecoder) .^ 2, 2) + 1;
% interference power [I_i^i] at the i-th private layer
intPowTerm = zeros(user);
% common stream power (1 layer)
comPowTerm = abs(bcChannel' * comPrecoder) .^ 2;
% total receive power [T^c] (i.e. total power at the first layer)
powTerm = priPowTerm(:, 1) + comPowTerm;

for iUser = 1 : user
    for iLayer = 1 : user
        if iLayer ~= 1
            % remaining private power at the i-th layer (i > 1)
            priPowTerm(iUser, iLayer) = priPowTerm(iUser, iLayer - 1) - abs(bcChannel(:, iUser)' * priPrecoder(:, iLayer - 1)) .^ 2;
        end
        % interference power at the i-th layer (when decoding common stream)
        intPowTerm(iUser, iLayer) = priPowTerm(iUser, iLayer) - abs(bcChannel(:, iUser)' * priPrecoder(:, iLayer)) .^ 2;
    end
end

% SIC stops after users decoding own layer
priPowTerm = tril(priPowTerm);
priPowTerm(priPowTerm == 0) = NaN;
intPowTerm = tril(intPowTerm);
intPowTerm(intPowTerm == 0) = NaN;

% optimum MMSE equalizers [g]
comEqualizer = zeros(user, 1);
priEqualizer = zeros(user);
for iUser = 1 : user
    comEqualizer(iUser) = comPrecoder' * bcChannel(:, iUser) / powTerm(iUser);
    for iLayer = 1 : user
        priEqualizer(iUser, iLayer) = priPrecoder(:, iLayer)' * bcChannel(:, iUser) / priPowTerm(iUser, iLayer);
    end
end

% MMSEs
comMmse = powTerm .\ priPowTerm(:, 1);
priMmse = priPowTerm .\ intPowTerm;

% optimum MMSE weights
comWeight = 1 ./ comMmse;
priWeight = 1 ./ priMmse;

% common rate must be achievable for all users (treat private part as noise)
comRate = min(real(log2(powTerm ./ priPowTerm(:, 1))));
% private rate (no influence from common part; rate of each layer must be achievable for those decode it)
priRate = min(real(log2(priPowTerm ./ intPowTerm)));

end
