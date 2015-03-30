% Function to fit the Hawkes parameters with MLE
% Author: Trey Huffine

function logLikelihood = HawkesMLE(parameters,times)
% Function used to in optmization to fit the Hawkes parameters.  Use the
% negative of the minimization to find the maximum.
% Input: 
%   parameters - The intial guess for mu, alpha, and beta being fitted
%   times - A list of trade arrival times
% Output:
%   logLikelihood - LL function used for fitting

% Equations:
% LL = -T(end)*mu + alpha/beta * sum(exp(-beta*(t(end)-ti)-1) + ...
%      sum(log(mu+alpha*R(i))
% R(i) = exp(-beta*(t_(i) - t_(i-1))) * (1+R(i-1))
% R(1) = 0
% firstSum = alpha/beta * sum(exp(-beta*(t(end)-ti)-1)
% secondSum = sum(log(mu+alpha*R(i))

mu = parameters(1);
alpha = parameters(2);
beta = parameters(3);
T = size(times,1);

% find firstSum
timeDifference = times(T) - times;
timeExponential = exp(-beta*timeDifference)-1;
firstSum = alpha/beta * sum(timeExponential);

% find secondSum
R = zeros(T,1);
for timeCounter = 2:T
    R(timeCounter) = exp(-beta * (times(timeCounter) - ...
        times(timeCounter-1)))*(1+R(timeCounter-1));
end
secondSum = sum( log(mu + alpha*R) );
    
logLikelihood = -(-mu*times(T) + firstSum + secondSum);