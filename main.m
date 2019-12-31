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
% load('data/data.mat');
figure('Name', sprintf('Achievable rate region comparison for %d-user %d-tx deployment, with \\gamma = %d and SNR = %d dB', user, tx, channelRelativeStrength, pow2db(snr)), 'Position', [500 200 900 600]);
legendString = cell(nAngles, 1);
for iAngle = 1 : nAngles
    subplot(2, 2, iAngle);
    % DPC
    dpcInstance = cell2mat(dpcRate(iAngle, :)');
    dpcInstance = remove_vertices(dpcInstance);
    dpcPlot = plot(dpcInstance(:, 1), dpcInstance(:, 2));
    legendString{1} = sprintf('DPC');
    hold on;

    % MU-LP
    mulpInstance = cell2mat(mulpRate(iAngle, :)');
    mulpInstance = remove_vertices(mulpInstance);
    mulpPlot = plot(mulpInstance(:, 1), mulpInstance(:, 2));
    legendString{2} = sprintf('MU-LP');
    hold on;

    % NOMA
    nomaInstance = cell2mat(nomaRate(iAngle, :)');
    nomaInstance = remove_vertices(nomaInstance);
    nomaPlot = plot(nomaInstance(:, 1), nomaInstance(:, 2));
    legendString{3} = sprintf('NOMA');
    hold on;

    % Single-layer RS
    slrsInstance = cell2mat(slrsRate(iAngle, :)');
    slrsInstance = remove_vertices(slrsInstance);
    slrsPlot = plot(slrsInstance(:, 1), slrsInstance(:, 2));
    legendString{4} = sprintf('SLRS');
    hold on;

    % RS
    rsInstance = cell2mat(rsRate(iAngle, :)');
    rsInstance = remove_vertices(rsInstance);
    rsPlot = plot(rsInstance(:, 1), rsInstance(:, 2));
    legendString{5} = sprintf('RS');

    hold off;
    grid on; grid minor;
    title(sprintf('\\theta = %d\\pi / 9', iAngle));
    legend(legendString, 'location', 'sw');
    xlabel('R_1 [bps/Hz]');
    ylabel('R_2 [bps/Hz]');
end
save([pwd '/data/data.mat']);
Push.pushNote(Push.Devices, 'MATLAB Assist', 'Job''s done!');
