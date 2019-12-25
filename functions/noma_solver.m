function [precoder, wsr] = noma_solver(weight, bcChannel, snr, equalizer, mmseWeight, order)
%NOMA_SOLVER Summary of this function goes here
%   Detailed explanation goes here


[tx, user] = size(bcChannel);

cvx_begin quiet

    variable precoder(tx, user) complex;

    % receive power terms [T]
    powTerm = cvx(zeros(user));

    % total receive power [T^c] (i.e. power at the first layer)
    powTerm(:, 1) = sum_square_abs(bcChannel' * precoder, 2) + 1;
    for iUser = 1 : user
        for iLayer = 2 : user
            % remaining power at the i-th layer (write as summation for cvx implementation)
            powTerm(iUser, iLayer) = sum_square_abs(bcChannel(:, iUser)' * precoder(:, order(iLayer : end)), 2) + 1;
        end
    end

    % MSEs and augmented WMSEs
    mse = cvx(zeros(user));
    wmse = cvx(zeros(user));
    % user rates on streams
    rate = cvx(zeros(user));
    for iUser = 1 : user
        for iLayer = 1 : user
            if isnan(equalizer(iUser, iLayer))
                % invalid stream (decoding already finished in previous layers)
                mse(iUser, iLayer) = NaN;
                wmse(iUser, iLayer) = NaN;
                rate(iUser, iLayer) = NaN;
            else
                mse(iUser, iLayer) = square_abs(equalizer(iUser, iLayer)) * powTerm(iUser, iLayer) - 2 * real(equalizer(iUser, iLayer) * bcChannel(:, iUser)' * precoder(:, iLayer)) + 1;
                wmse(iUser, iLayer) = mmseWeight(iUser, iLayer) * mse(iUser, iLayer) - log2(mmseWeight(iUser, iLayer));
                rate(iUser, iLayer) = 1 - wmse(iUser, iLayer);
            end
        end
    end

    % calculate weighted sum-rate
    wsr = 0;
    for iUser = 1 : user
        % layer rates for a specific user
        layerRate = rate(iUser, :);
        % CVX does not support comparison including NaN, so we remove invalid entries before minimization
        clsIdx = cvx_classify(layerRate);
        % 13 -> invalid
        wsr = wsr + weight(iUser) * min(layerRate(clsIdx ~= 13));
    end

    % solve weighted sum-rate maximization problem
    maximize wsr;
    subject to
        precoder(:)' * precoder(:) <= snr;

cvx_end

end
