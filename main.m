initialize; config;
Push.pushNote(Push.Devices, 'MATLAB Assist', sprintf('''%s'' is running', mfilename));

%% R-E region calculation
nAngles = size(channelRelativeAngle, 2);
nWeights = size(weight, 2);

dpcRate = cell(nAngles, nWeights);
mulpRate = cell(nAngles, nWeights);
nomaRate = cell(nAngles, nWeights);
slrsRate = cell(nAngles, nWeights);
rsRate = cell(nAngles, nWeights);

try
    for iAngle = 1 : nAngles
        % update BC channel of user 2
        bcChannel(:, :, 2) = kron(channelRelativeStrength, exp(1j * channelRelativeAngle(iAngle) * (0 : 3)));
        for iWeight = 1 : nWeights
            [dpcRate{iAngle, iWeight}] = dpc_rate(weight(:, iWeight), bcChannel, snr, tolerance);
            [mulpRate{iAngle, iWeight}] = mulp_rate(weight(:, iWeight), bcChannel, snr, tolerance);
            [nomaRate{iAngle, iWeight}] = noma_rate(weight(:, iWeight), bcChannel, snr, tolerance);
            [slrsRate{iAngle, iWeight}] = slrs_rate(weight(:, iWeight), bcChannel, snr, tolerance, rsRatio);
            [rsRate{iAngle, iWeight}] = rs_rate(weight(:, iWeight), bcChannel, snr, tolerance, rsRatio);
        end
    end
catch
    Push.pushNote(Push.Devices, 'MATLAB Assist', 'Houston, we have a problem');
end
%% R-E region comparison
figure('Name', sprintf('Achievable rate region comparison for %d-user %d-tx deployment, with \gamma = %d and SNR = %d dB', [user, rx, channelRelativeStrength, snr]));
legendString = cell(nAngles, 1);
for iAngle = 1 : nAngles
    subplot(2, 2, iAngle);
    % DPC
    dpcPlot = plot(cell2mat(dpcRate(iAngle, :)'));
    legendString{1} = sprintf('DPC');
    hold on;

    % MU-LP
    mulpPlot = plot(cell2mat(mulpRate(iAngle, :)'));
    legendString{2} = sprintf('MU-LP');
    hold on;

    % NOMA
    nomaIdx = convhull(cell2mat(nomaRate(iAngle, :)'));
    nomaPlot = plot(sortrow(cell2mat(nomaRate(iAngle, nomaIdx)')));
    legendString{3} = sprintf('NOMA');
    hold on;

    % Single-layer RS
    slrsPlot = plot(cell2mat(slrsRate(iAngle, :)'));
    legendString{4} = sprintf('SLRS');
    hold on;

    % RS
    rsIdx = convhull(cell2mat(rspRate(iAngle, :)'));
    rsPlot = plot(sortrow(cell2mat(rsRate(iAngle, rsIdx)')));
    legendString{5} = sprintf('RS');

    hold off;
    grid on; grid minor;
    legend(legendString);
    xlabel('R_1 [bps/Hz]');
    ylabel('R_2 [bps/Hz]');
    save([pwd '/data/data.mat']);
end
Push.pushNote(Push.Devices, 'MATLAB Assist', 'Job''s done!');
