function [rate] = mulp_rate(weight, bcChannel, snr, tolerance)
%MULP_RATE Summary of this function goes here
%   Detailed explanation goes here


[rx, tx, user] = size(bcChannel);

% initialize precoders [p] (tx * rx * user) by matched filters
precoder = sqrt(snr / user) * conj(bcChannel) ./ vecnorm(bcChannel);
precoder = permute(precoder, [2, 1, 3]);

% noise covariance matrix [R_vv] (rx * rx * user)
noiseCovMat = repmat(eye(rx), [1, 1, user]);
% MMSE receive filter [A] (user * rx)
mmseEqualizer = zeros(user, rx);
% MSE covariance matrix [E_k] (rx * rx * user)
mseCovMat = zeros(rx, rx, user);
% user rates
rate = zeros(1, user);

for iUser = 1:user
    for jUser = 1:user
        if jUser ~= iUser
            noiseCovMat(:, :, iUser) = noiseCovMat(:, :, iUser) + bcChannel(:, :, iUser) * precoder(:, :, jUser) * precoder(:, :, jUser)' * bcChannel(:, :, iUser)';
        end
    end
    mmseEqualizer(iUser, :) = precoder(:, :, iUser)' * bcChannel(:, :, iUser)' / (bcChannel(:, :, iUser) * precoder(:, :, iUser) * precoder(:, :, iUser)' * bcChannel(:, :, iUser)' + noiseCovMat(:, :, iUser));
    mseCovMat(:, :, iUser) = inv(eye(rx) + precoder(:, :, iUser)' * bcChannel(:, :, iUser)' / noiseCovMat(:, :, iUser) * bcChannel(:, :, iUser) * precoder(:, :, iUser));
    rate(iUser) = log2(det(mseCovMat(:, :, iUser)));
end

end
