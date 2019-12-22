initialize; config;

nAngles = size(channelRelativeAngle, 2);
nWeights = size(weight, 2);

dpcRate = cell(nAngles, nWeights);
rsRate = cell(nAngles, nWeights);
mulpRate = cell(nAngles, nWeights);
for iAngle = 1:nAngles
    % update BC channel of user 2
    bcChannel(:, :, 2) = kron(channelRelativeStrength, exp(1j * channelRelativeAngle(iAngle) * (0:3)));
    for iWeight = 1:nWeights
%         [dpcRate{iAngle, iWeight}] = dpc_rate(weight(:, iWeight), bcChannel, snr, tolerance);
%         [rsRate{iAngle, iWeight}] = rs_rate(weight(:, iWeight), bcChannel, snr, tolerance, rsRatio);
        [mulpRate] = mulp_rate(weight(:, iWeight), bcChannel, snr, tolerance);
    end
end
