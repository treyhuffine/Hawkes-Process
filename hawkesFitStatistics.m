function [ totalRMSE, AIC, BIC ] = hawkesFitStatistics( empirical, hawkesModel, logLikelihood )
%Metrics to quantify the fit of MLE
p = 3; % mu, alpha, beta
n = numel(empirical);

[totalRMSE, ~] = RMSE(empirical, hawkesModel); 

% AIC and BIC are normalized by number of data points in window
AIC = -2*logLikelihood/n + 2*p/n;
BIC = -2*logLikelihood/n + log(n)*p/n;

end

