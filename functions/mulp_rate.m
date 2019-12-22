function [rate] = mulp_rate(weight, bcChannel, snr, tolerance)
% Function:
%   - compute the achievable user rates with multi-user linear precoding
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
%   - maximize weighted sum-rate using weighted MMSE
%   - choose MMSE weights to ensure the WMMSE-gradient is identical to the WSR-gradient
%   - alternate between WMMSE precoder optimization and MMSE weight update: (equalizer, MMSE weights) | precoder -> precoder | (equalizer, MMSE weights)
%
% Reference(s):
%   - S. S. Christensen, R. Agarwal, E. D. Carvalho, and J. M. Cioffi, "Weighted sum-rate maximization using weighted MMSE for MIMO-BC beamforming design,"”" IEEE Transactions on Wireless Communications, vol. 7, no. 12, pp. 4792–4799, Dec 2008.
%
% Author & Date: Yang (i@snowztail.com) - 22 Dec 19


[rx, tx, user] = size(bcChannel);

% row-wise stacked channel matrix [H] ((rx * user) * tx)
stackedChannel = reshape(permute(bcChannel, [1, 3, 2]), [user * rx, tx]);

% initialize precoders [B] (tx * rx * user) by matched filter
wmmsePrecoder = sqrt(snr / user) * conj(bcChannel) ./ vecnorm(bcChannel);
wmmsePrecoder = permute(wmmsePrecoder, [2, 1, 3]);

% initialize noise and MSE covariance matrices
[noiseCovMat, mmseCovMat, ~] = covariance(bcChannel, wmmsePrecoder);

isConverged = false;
capacity = 0;
while (~isConverged)
    % MMSE equalizer [A] {1 * user}
    mmseEqualizer = cell(1, user);
    % MMSE weight [W] {1 * user}
    mmseWeight = cell(1, user);

    % update MMSE equalizers and weights based on previous precoders
    for iUser = 1:user
        mmseEqualizer{iUser} = wmmsePrecoder(:, :, iUser)' * bcChannel(:, :, iUser)' / (bcChannel(:, :, iUser) * wmmsePrecoder(:, :, iUser) * wmmsePrecoder(:, :, iUser)' * bcChannel(:, :, iUser)' + noiseCovMat(:, :, iUser));
        % select MMSE weight to ensure WMMSE-gradient equals WSR-gradient
        mmseWeight{iUser} = weight(iUser) / mmseCovMat(:, :, iUser);
    end

    % construct block-diagonal matrices
    mmseWeight = blkdiag(mmseWeight{:});
    mmseEqualizer = blkdiag(mmseEqualizer{:});

    % update WMMSE precoders [B] (tx * (rx * user))
    wmmsePrecoder = (stackedChannel' * mmseEqualizer' * mmseWeight * mmseEqualizer * stackedChannel + trace(mmseWeight * (mmseEqualizer * mmseEqualizer')) * eye(tx) / snr) \ stackedChannel' * mmseEqualizer' * mmseWeight;
    % adapt precoder to transmit power constraint
    wmmsePrecoder = sqrt(snr / trace(wmmsePrecoder * wmmsePrecoder')) * wmmsePrecoder;
    % reshape precoder as (tx * rx * user)
    wmmsePrecoder = reshape(wmmsePrecoder, [tx, rx, user]);
    % update covariance matrices
    [noiseCovMat, mmseCovMat, rate] = covariance(bcChannel, wmmsePrecoder);

    if (sum(rate) - capacity) / sum(rate) <= tolerance
        isConverged = true;
    end
    capacity = sum(rate);
end

end


function [noiseCovMat, mmseCovMat, rate] = covariance(bcChannel, wmmsePrecoder)
    [rx, ~, user] = size(bcChannel);
    % noise covariance matrix [R_vv] (rx * rx * user)
    noiseCovMat = repmat(eye(rx), [1, 1, user]);
    % MSE covariance matrix [E_k] (rx * rx * user)
    mmseCovMat = zeros(rx, rx, user);
    % user rate (1 * user)
    rate = zeros(1, user);
    for iUser = 1:user
        for jUser = 1:user
            if jUser ~= iUser
                noiseCovMat(:, :, iUser) = noiseCovMat(:, :, iUser) + bcChannel(:, :, iUser) * wmmsePrecoder(:, :, jUser) * wmmsePrecoder(:, :, jUser)' * bcChannel(:, :, iUser)';
            end
        end
        mmseCovMat(:, :, iUser) = inv(eye(rx) + wmmsePrecoder(:, :, iUser)' * bcChannel(:, :, iUser)' / noiseCovMat(:, :, iUser) * bcChannel(:, :, iUser) * wmmsePrecoder(:, :, iUser));
        rate(iUser) = real(log2(det(inv(mmseCovMat(:, :, iUser)))));
    end
end
