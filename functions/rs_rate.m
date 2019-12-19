function [rate] = rs_rate(weight, bcChannel, snr, tolerance, rsRatio)
%RS_RATE Summary of this function goes here
%   Detailed explanation goes here
% MISO

[rx, tx, user] = size(bcChannel);
% reshape BC channel matrix [H] (tx * rx * user) with Ref 1
% bcChannel = conj(permute(bcChannel, [2, 1, 3]));
bcChannel = squeeze(permute(bcChannel, [2, 1, 3]));

% initialize common precoder using MRT & SVD (see Ref 2)
[u, ~, ~] = svd(bcChannel);
% largest left singular vector of channel matrix
llsv = u(:, 1);
% common precoder (tx * 1)
comPrecoder = sqrt(snr * (1 - rsRatio)) * llsv;
clearvars u;
% private precoder (tx * user)
priPrecoder = sqrt(snr * (rsRatio / user)) * bcChannel ./ vecnorm(bcChannel);

% private stream power [T]
powPriTerm = sum(abs(bcChannel' * priPrecoder) .^ 2, 2) + 1;
% common stream power
powComTerm = abs(bcChannel' * comPrecoder) .^ 2;
% total receive power (common + private)
powTerm = powPriTerm + powComTerm;
% interference power
powIntTerm = zeros(user, 1);
for iUser = 1:user
    powIntTerm(iUser) = powPriTerm(iUser) - abs(bcChannel(:, iUser)' * priPrecoder(:, iUser)) .^ 2;
end

end
