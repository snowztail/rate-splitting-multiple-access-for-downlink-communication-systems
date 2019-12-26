function [rate] = dpc_rate(weight, bcChannel, snr, tolerance)
% Function:
%   - compute the achievable user rates with dirty paper coding
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
%   - achieves the channel capacity, no power penalty
%   - the receiver behaves as if no interference exists
%   - one multiuser rate point at a time
%   - encoding order depends on user weight
%   - obtain rate region by varying user weights (we can also examine all possible encoding orders)
%
% Reference(s):
%   - H. Viswanathan, S. Venkatesan, and H. Huang, "Downlink capacity evaluation of cellular networks with known-interference cancellation," IEEE Journal on Selected Areas in Communications, vol. 21, no. 5, pp. 802-811, June 2003.
%
% Author & Date: Yang (i@snowztail.com) - 18 Dec 19

[rx, tx, user] = size(bcChannel);
% dual MAC channel
macChannel = conj(permute(bcChannel, [2, 1, 3]));
% uplink covariance matrix [Q]
covMat = zeros(rx, rx, user);
% determine serving priority by weight
[weight, priority] = sort(weight, 'descend');
% sort channels
bcChannel = bcChannel(:, :, priority);
macChannel = macChannel(:, :, priority);

% padding weight for simple implementation
weight(end + 1) = 0;

isConverged = false;
capacity = 0;
while (~isConverged)
    % gradient of the objective function
    gradient = zeros(rx, rx, user);
    for iUser = 1 : user
        for iSucUser = iUser : user
            gradient(:, :, iUser) = gradient(:, :, iUser) + (weight(iSucUser) - weight(iSucUser + 1)) *...
            (bcChannel(:, :, iUser) / (eye(tx) + sum(macChannel(:, :, 1 : iSucUser) .* covMat(:, :, 1 : iSucUser) .* bcChannel(:, :, 1 : iSucUser), 3)) * macChannel(:, :, iUser));
        end
    end

    % principal eigenvalue and eigenvectors of gradient
    eigVal = zeros(1, user);
    eigVec = zeros(rx, user);
    for iUser = 1:user
        [V, D] = eig(gradient(:, :, iUser));
        [eigVal(iUser), idx] = max(diag(D));
        eigVec(:, iUser) = V(:, idx);
    end
    clearvars V D idx;

    % obtain optimal user index
    [~, usrIdx] = max(eigVal);

    % objective function with variable [t]
    objVal = @(ratio) -objective(ratio, covMat, usrIdx, eigVec, snr, weight, bcChannel, macChannel);
    % obtain optimal solution
    ratioOpt = fminbnd(objVal, 0, 1);

    % update covariance matrix
    for iUser = 1:user
        covMat(:, :, iUser) = ratioOpt * covMat(:, :, iUser);
    end
    covMat(:, :, usrIdx) = covMat(:, :, usrIdx) + (1 - ratioOpt) * snr * eigVec(:, usrIdx) * eigVec(:, usrIdx)';

    % calculate user rates
    rate = zeros(1, user);
    for iUser = 1 : user
        if iUser == 1
            rate(priority(iUser)) = log2(det(eye(tx) + sum(macChannel(:, :, 1 : iUser) .* covMat(:, :, 1 : iUser) .* bcChannel(:, :, 1 : iUser), 3)));
        else
            rate(priority(iUser)) = log2(det(eye(tx) + sum(macChannel(:, :, 1 : iUser) .* covMat(:, :, 1 : iUser) .* bcChannel(:, :, 1 : iUser), 3)) / det(eye(tx) + sum(macChannel(:, :, 1 : iUser - 1) .* covMat(:, :, 1 : iUser - 1) .* bcChannel(:, :, 1 : iUser - 1), 3)));
        end
    end
    if (sum(rate) - capacity) / sum(rate) <= tolerance
        isConverged = true;
    end
    capacity = sum(rate);
end

end


%% Subfunctions
% objective function [f]
function [objVal] = objective(ratio, covMat, usrIdx, eigVec, snr, weight, bcChannel, macChannel)
    [~, tx, user] = size(bcChannel);
    for iUser = 1 : user
        covMat(:, :, iUser) = ratio * covMat(:, :, iUser);
    end
    % optimal user
    covMat(:, :, usrIdx) = covMat(:, :, usrIdx) + (1 - ratio) * snr * eigVec(:, usrIdx) * eigVec(:, usrIdx)';
    % objective function
    objVal = 0;
    for iUser = 1 : user
        for iPreUser = 1 : iUser
            objVal = objVal + (weight(iPreUser) - weight(iPreUser + 1)) * log2(det(eye(tx) + sum(macChannel(:, :, 1 : iPreUser) .* covMat(:, :, 1 : iPreUser) .* bcChannel(:, :, 1 : iPreUser), 3)));
        end
    end
end
