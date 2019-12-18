function [rate] = dpc_rate(weight, bcChannel, snr)
%DPC_RATE Summary of this function goes here
%   Detailed explanation goes here

[rx, tx, user] = size(bcChannel);
% dual MAC channel
macChannel = conj(permute(bcChannel, [2, 1, 3]));
% uplink covariance matrix [\Q]
covMat = zeros(rx, rx, user);
% gradient of the objective function
weight(end + 1) = 0;
for iUser = 1:user
    for iSuccUser = iUser:user
        gradFun = gradFun + (weight(iSuccUser) - weight(iSuccUser + 1)) * (bcChannel(:, :, iUser) / (eye(tx) + sum(macChannel(:, :, 1:iSuccUser) * covMat(:, :, 1:iSuccUser) * bcChannel(:, :, 1:iSuccUser))) * macChannel(:, :, iUser));
    end
end
weight = weight(1: end - 1);

end
