function [comPrecoder, priPrecoder, wsr] = rs_solver(weight, bcChannel, snr, comEqualizer, priEqualizer, comWeight, priWeight)
% Function:
%   - solve the optimum common and private precoders for 1-layer rate-splitting multiple access
%
% InputArg(s):
%   - weight [u] (user * 1): user weights
%   - bcChannel [H] (rx * tx * user): broadcast channel response
%   - snr [\rho]: signal-to-noise ratio, which equals transmit power since unit noise power assumed
%   - comEqualizer [g_i^c] (user * 1): optimum MMSE equalizers for common stream
%   - priEqualizer [g_i^i] (user * 1): optimum MMSE equalizers for (self) private stream
%   - comWeight [u_i^c] (user * 1): optimum MMSE weights for common stream
%   - priWeight [u_i^i] (user * 1): optimum MMSE weights for (self) private stream
%
% OutputArg(s):
%   - comPrecoder [p_i^c] (tx * 1): optimum common precoder
%   - priPrecoder [p_i^i] (tx * user): optimum private precoders
%   - wsr: achievable weighted sum rate
%
% Comment(s):
%   - for 1-layer RS on MU-MISO systems only
%   - assume common rate is with unit weight
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 20 Dec 19

% for MU-MISO channels
[tx, user] = size(bcChannel);

cvx_begin quiet

    variable comPrecoder(tx, 1) complex;
    variable priPrecoder(tx, user) complex;

    % private stream power [T_i^i, I_i^c]
    priPowTerm = sum_square_abs(bcChannel' * priPrecoder, 2) + 1;
    % total receive power [T_i^c] (common + private)
    powTerm = square_abs(bcChannel' * comPrecoder) + priPowTerm;

    % MSEs
    comMse = cvx(user, 1);
    priMse = cvx(user, 1);
    for iUser = 1 : user
        comMse(iUser) = square_abs(comEqualizer(iUser)) * powTerm (iUser) - 2 * real(comEqualizer(iUser) * bcChannel(:, iUser)' * comPrecoder) + 1;
        priMse(iUser) = square_abs(priEqualizer(iUser)) * priPowTerm(iUser) - 2 * real(priEqualizer(iUser) * bcChannel(:, iUser)' * priPrecoder(:, iUser)) + 1;
    end

    % augmented WMSEs
    comWmse = comWeight .* comMse - log2(comWeight);
    priWmse = priWeight .* priMse - log2(priWeight);

    % total user rate (common + private)
    comRate = 1 - comWmse;
    priRate = 1 - priWmse;
    % weighted sum-rate (assume common rate is with unit weight)
    wsr = min(comRate) + sum(weight .* priRate);

    % solve weighted sum-rate maximization problem
    maximize wsr;
    subject to
        priPrecoder(:)' * priPrecoder(:) + comPrecoder' * comPrecoder <= snr;

cvx_end

end
