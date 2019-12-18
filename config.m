%% Transceiver
% number of transmit/receive antennas
tx = 4;
rx = 1;
%% SNR
snr = db2pow(20);

%% User
% number of users
user = 2;
% weight (user * instance)
weight(2, :) = exp([-3, -1:0.05:1, 3]);
weight(1, :) = 1;

%% Channel
% angle between the channels of user 1 and 2 [\theta]
channelRelativeAngle = (1:4) * pi / 9;
% channel strength ratio [\gamma]
channelRelativeStrength = 1;
% broadcast channel gains (rx * tx * user)
bcChannel(:, :, 2) = kron(channelRelativeStrength, exp(1j * (0:3)));
bcChannel(:, :, 1) = ones(rx, tx);

[rate] = dpc_rate(weight, bcChannel, snr);
