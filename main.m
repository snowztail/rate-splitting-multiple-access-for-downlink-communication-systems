initialize; config;

for iAngle = 1:size(channelRelativeAngle)
    % update BC channel of user 2
    bcChannel(:, :, 2) = kron(channelRelativeStrength, exp(1j * channelRelativeAngle(iAngle) * (0:3)));
    % [rate] = dpc_rate(weight(:, 1), bcChannel, snr, tolerance);
    [rate] = rs_rate(weight(:, 1), bcChannel, snr, tolerance, rsRatio);
end
