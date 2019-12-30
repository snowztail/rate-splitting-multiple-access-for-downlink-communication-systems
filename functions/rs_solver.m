function [comPrecoder, priPrecoder, wsr] = rs_solver(weight, bcChannel, snr, comEqualizer, priEqualizer, comWeight, priWeight)
% Function:
%   - solve the optimum common and private precoders for general rate-splitting multiple access
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
%   - row -> user, column -> stream; (i, j) entry means the j-th stream decoded by user-i (NaN if invalid)
%   - require sorted terms based on decoding order
%   - user-i first decode the common stream, then decode the private streams 1 to i (self stream)
%   - the common rate goes to user-1 since the common channel is aligned to user-1
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 27 Dec 19


[tx, user] = size(bcChannel);

cvx_begin quiet

    variable comPrecoder(tx, 1) complex;
    variable priPrecoder(tx, user) complex;

    % private stream power [T_i^i, I_i^c] ($user$ layers)
    priPowTerm = cvx(zeros(user));
    % total private power at the first layer
    priPowTerm(:, 1) = sum_square_abs(bcChannel' * priPrecoder, 2) + 1;
    % common stream power (1 layer)
    comPowTerm = square_abs(bcChannel' * comPrecoder);
    % total receive power [T^c] (i.e. total power at the first layer)
    powTerm = priPowTerm(:, 1) + comPowTerm;

    for iUser = 1 : user
        for iLayer = 2 : user
            % remaining private power at the i-th layer (i > 1)
            priPowTerm(iUser, iLayer) = sum_square_abs(bcChannel(:, iUser)' * priPrecoder(:, iLayer : end), 2) + 1;
        end
    end

    % MSEs
    comMse = cvx(zeros(user, 1));
    priMse = cvx(zeros(user));
    % augmented WMSEs
    comWmse = cvx(zeros(user, 1));
    priWmse = cvx(zeros(user));
    % user rates
    comRate = cvx(zeros(user, 1));
    priRate = cvx(zeros(user));
    for iUser = 1 : user
        % common stream: MSE -> WMSE -> rate
        comMse(iUser) = square_abs(comEqualizer(iUser)) * powTerm (iUser) - 2 * real(comEqualizer(iUser) * bcChannel(:, iUser)' * comPrecoder) + 1;
        comWmse(iUser) = comWeight(iUser) * comMse(iUser) - log2(comWeight(iUser));
        comRate(iUser) = 1 - comWmse(iUser);
        for iLayer = 1 : user
            if isnan(priEqualizer(iUser, iLayer))
                % invalid stream (decoding already finished in previous layers)
                priMse(iUser, iLayer) = NaN;
                priWmse(iUser, iLayer) = NaN;
                priRate(iUser, iLayer) = NaN;
            else
                priMse(iUser, iLayer) = square_abs(priEqualizer(iUser, iLayer)) * priPowTerm(iUser, iLayer) - 2 * real(priEqualizer(iUser, iLayer) * bcChannel(:, iUser)' * priPrecoder(:, iLayer)) + 1;
                priWmse(iUser, iLayer) = priWeight(iUser, iLayer) * priMse(iUser, iLayer) - log2(priWeight(iUser, iLayer));
                priRate(iUser, iLayer) = 1 - priWmse(iUser, iLayer);
            end
        end
    end

    % common rate should be achievable for all users (goes to user-1 as the common channel is aligned to user-1)
    comWsr = weight(1) * min(comRate);
    % weighted-sum of private rates is obtained by <1> select layer rates which can be achieved by those need to decode it (by min over valid entries) <2> find the corresponding user weight by decoding order
    priWsr = 0;
    for iLayer = 1 : user
        % user private rates for a specific layer
        layerRate = priRate(:, iLayer);
        % CVX does not support comparison including NaN, so we remove invalid entries before minimization
        clsIdx = cvx_classify(layerRate);
        % 13 -> invalid
        priWsr = priWsr + weight(iLayer) * min(layerRate(clsIdx ~= 13));
    end
    % total weighted sum-rate
    wsr = comWsr + priWsr;

    % solve weighted sum-rate maximization problem
    maximize wsr;
    subject to
        priPrecoder(:)' * priPrecoder(:) + comPrecoder' * comPrecoder <= snr;

cvx_end

end
