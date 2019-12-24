function [comEqualizer, priEqualizer, comWeight, priWeight, comRate, priRate] = rs_terms(bcChannel, comPrecoder, priPrecoder)
% Function:
%   - construct power terms, equalizers, MMSEs, and weights for precoder optimization
%   - calculate achievable common and private rates
%
% InputArg(s):
%   - bcChannel [H] (rx * tx * user): broadcast channel response
%   - comPrecoder [p_i^c] (tx * 1): optimum common precoder
%   - priPrecoder [p_i^i] (tx * user): optimum private precoders
%   - snr [\rho]: signal-to-noise ratio, which equals transmit power since unit noise power assumed
%
% OutputArg(s):
%   - comEqualizer [g_i^c] (user * 1): optimum MMSE equalizers for common stream
%   - priEqualizer [g_i^i] (user * 1): optimum MMSE equalizers for (self) private stream
%   - comWeight [u_i^c] (user * 1): optimum MMSE weights for common stream
%   - priWeight [u_i^i] (user * 1): optimum MMSE weights for (self) private stream
%   - comRate [R_c]: achievable common rate
%   - priRate [R_i] (user * 1): user private rates
%
% Comment(s):
%   - for 1-layer RS on MU-MISO systems only
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 20 Dec 19


[~, user] = size(bcChannel);
% private stream power [T_i^i, I_i^c]
priPowTerm = sum(abs(bcChannel' * priPrecoder) .^ 2, 2) + 1;
% total receive power [T_i^c] (common + private)
powTerm = abs(bcChannel' * comPrecoder) .^ 2 + priPowTerm;
% interference power [I_i^i]
intPowTerm = zeros(user, 1);
for iUser = 1 : user
    intPowTerm(iUser) = priPowTerm(iUser) - abs(bcChannel(:, iUser)' * priPrecoder(:, iUser)) .^ 2;
end

% optimum MMSE equalizers [g]
comEqualizer = zeros(user, 1);
priEqualizer = zeros(user, 1);
for iUser = 1 : user
    comEqualizer(iUser) = comPrecoder' * bcChannel(:, iUser) / powTerm(iUser);
    priEqualizer(iUser) = priPrecoder(:, iUser)' * bcChannel(:, iUser) / priPowTerm(iUser);
end

% MMSEs
comMmse = powTerm .\ priPowTerm;
priMmse = priPowTerm .\ intPowTerm;

% optimum MMSE weights
comWeight = 1 ./ comMmse;
priWeight = 1 ./ priMmse;

% common rate must be achievable for all users (treat private part as noise)
comRate = min(log2(powTerm ./ priPowTerm));
% private rate (no influence from common part)
priRate = log2(priPowTerm ./ intPowTerm);

end
