function [precoder, wsr] = noma_solver(weight, bcChannel, snr, equalizer, mmseWeight)
%NOMA_SOLVER Summary of this function goes here
%   Detailed explanation goes here


[tx, user] = size(bcChannel);

cvx_begin quiet

    variable precoder(tx, user) complex;

    % receive power terms [T]
    powTerm = cvx(user);
    % interference power terms [I]
    intPowTerm = cvx(user);

    % total receive power [T^c] (i.e. power at the first layer)
    powTerm(:, 1) = sum(abs(bcChannel' * precoder) .^ 2, 2) + 1;
    powTerm(:, 1) = sum_square_abs(bcChannel' * precoder, 2) + 1;
    for iUser = 1 : user
        for iLayer = 1 : user
            if iLayer == 1
                % interference power at the first layer
                intPowTerm(iUser, iLayer) = powTerm(iUser, iLayer) - square_abs(bcChannel(:, iUser)' * precoder(:, order(iLayer)));
            else
                % remaining power at the i-th layer
                powTerm(iUser, iLayer) = powTerm(iUser, iLayer - 1) - square_abs(bcChannel(:, iUser)' * precoder(:, order(iLayer - 1)));
                % interference power at the i-th layer
                intPowTerm(iUser, iLayer) = intPowTerm(iUser, iLayer - 1) - square_abs(bcChannel(:, iUser)' * precoder(:, order(iLayer)));
            end
        end
    end

    % MSEs
    mse = cvx(user);
    for iUser = 1 : user
        for iLayer = 1 : user
            mse(iUser, iLayer) = square_abs(equalizer(iUser, iLayer)) * powTerm(iUser, iLayer) - 2 * real(equalizer(iUser, iLayer) * bcChannel(:, iUser)' * priPrecoder(:, order(iLayer - 1))) + 1;
        end
        % comMse(iUser) = square_abs(comEqualizer(iUser)) * powTerm (iUser) - 2 * real(comEqualizer(iUser) * bcChannel(:, iUser)' * comPrecoder) + 1;
        % priMse(iUser) = square_abs(priEqualizer(iUser)) * priPowTerm(iUser) - 2 * real(priEqualizer(iUser) * bcChannel(:, iUser)' * priPrecoder(:, iUser)) + 1;
    end

    % augmented WMSEs
    comWmse = comWeight .* comMse - log2(comWeight);
    priWmse = priWeight .* priMse - log2(priWeight);

    % total user rate (common + private)
    comRate = 1 - comWmse;
    priRate = 1 - priWmse;
    % weighted sum-rate (assume common rate is with unit weight)
    wsr = min(comRate) + sum(weight .* priRate);

    % solve weighted sum-rate maximization problem
    maximize wsr;
    subject to
        priPrecoder(:)' * priPrecoder(:) + comPrecoder' * comPrecoder <= snr;

cvx_end
end
