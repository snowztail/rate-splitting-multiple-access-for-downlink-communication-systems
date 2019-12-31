%% Transceiver
% number of transmit/receive antennas
tx = 4;
rx = 1;

%% SNR
snr = db2pow(20);
%% tolerance for convergence (percentage)
tolerance = 1e-6;

%% User
% number of users
user = 2;
% weight (user * instance)
weight(2, :) = 10 .^ ([-3, -1 : 0.05 : 1, 3]);
weight(1, :) = 1;

%% Channel
% angle between the channels of user 1 and 2 [\theta]
channelRelativeAngle = (1 : 4) * pi / 9;
% channel strength ratio [\gamma]
channelRelativeStrength = 1;
% broadcast channel gains (rx * tx * user)
bcChannel(:, :, 1) = ones(rx, tx);

%% RSMA
% power ratio for private message streams [\alpha]
rsRatio = 0.2;

%% Pushbullet (for notifications)
apiKey = 'o.KLFeL1TWlJeMit7JOTmJKg7DTfsxvaXQ';
Push = Pushbullet(apiKey);
