function [totalRMSE, fullDiff] = RMSE(empirical, hawkesModel)
    % Take matrices with the returns you are interested in.  The columns of
    % the empirical0 will be compared with the fitted model
[L W] = size(empirical);

% Finds the difference (NOT rmse) between the at each point regardless
% of the data
fullDiff = empirical-hawkesModel;
fullDiff(isnan(fullDiff)) = 0;

%  Find the RMSE
totalRMSE = sqrt(mean(fullDiff.^2));

