function [precoder, wsr] = noma_solver(weight, bcChannel, snr, equalizer, mmseWeight)
% Function:
%   - solve the optimum precoder (regarding weighted-sum rate) for non-orthogonal multiple access with successive interference cancellation
%
% InputArg(s):
%   - weight [u] (user * 1): user weights
%   - bcChannel [H] (rx * tx * user): broadcast channel response
%   - snr [\rho]: signal-to-noise ratio, which equals transmit power since unit noise power assumed
%   - equalizer [g] (user * 1): optimum MMSE equalizer
%   - mmseWeight [u^mmse] (user * 1): optimum MMSE weights
%
% OutputArg(s):
%   - precoder [p] (tx * user): optimum precoders maximizing WSR
%   - wsr: achievable weighted sum rate
%
% Comment(s):
%   - row -> user, column -> stream; (i, j) entry means the j-th stream decoded by user-i (NaN if invalid)
%   - require sorted terms based on decoding order
%   - user-i decodes layer 1 to i (terminates decoding once reaching self stream i)
%
% Reference(s):
%   - Y. Mao, B. Clerckx, and V. O. Li, "Rate-splitting multiple access for downlink communication systems: bridging, generalizing, and outperforming SDMA and NOMA," EURASIP Journal on Wireless Communications and Networking, vol. 2018, no. 1, 2018.
%
% Author & Date: Yang (i@snowztail.com) - 25 Dec 19


[tx, user] = size(bcChannel);

cvx_begin quiet

    variable precoder(tx, user) complex;

    % receive power terms [T]
    powTerm = cvx(zeros(user));

    % total receive power [T^c] (i.e. power at the first layer)
    powTerm(:, 1) = sum_square_abs(bcChannel' * precoder, 2) + 1;
    for iUser = 1 : user
        for iLayer = 2 : user
            % remaining power at the i-th layer (write as summation for cvx implementation)
            powTerm(iUser, iLayer) = sum_square_abs(bcChannel(:, iUser)' * precoder(:, iLayer : end), 2) + 1;
        end
    end

    % MSEs and augmented WMSEs
    mse = cvx(zeros(user));
    wmse = cvx(zeros(user));
    % user rates on streams
    rate = cvx(zeros(user));
    for iUser = 1 : user
        for iLayer = 1 : user
            if isnan(equalizer(iUser, iLayer))
                % invalid stream (decoding already finished in previous layers)
                mse(iUser, iLayer) = NaN;
                wmse(iUser, iLayer) = NaN;
                rate(iUser, iLayer) = NaN;
            else
                mse(iUser, iLayer) = square_abs(equalizer(iUser, iLayer)) * powTerm(iUser, iLayer) - 2 * real(equalizer(iUser, iLayer) * bcChannel(:, iUser)' * precoder(:, iLayer)) + 1;
                wmse(iUser, iLayer) = mmseWeight(iUser, iLayer) * mse(iUser, iLayer) - log2(mmseWeight(iUser, iLayer));
                rate(iUser, iLayer) = 1 - wmse(iUser, iLayer);
            end
        end
    end

    % calculate weighted sum-rate
    wsr = 0;
    for iLayer = 1 : user
        % user rates for a specific layer
        layerRate = rate(:, iLayer);
        % CVX does not support comparison including NaN, so we remove invalid entries before minimization
        clsIdx = cvx_classify(layerRate);
        % 13 -> invalid
        wsr = wsr + weight(iLayer) * min(layerRate(clsIdx ~= 13));
    end

    % solve weighted sum-rate maximization problem
    maximize wsr;
    subject to
        precoder(:)' * precoder(:) <= snr;

cvx_end

end
