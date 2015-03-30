%------------------------------------------------------------------------
% Fit a Hawkes Process using MLE
% Author: Trey Huffine
% University of California, Berkeley
% Master of Financial Engineering
%------------------------------------------------------------------------

close all, clear all, clc;
warning off;
tic %used for timing the algorith

%% Declare User Input
numTrades = 10000; % sets the window length
saveImages = false; % set to true to save an image each iteration 
% set pathmame to save images to
pathname = 'C:\your_path';
file = 'file_name.csv';
all_trades = csvread(file);
times = all_trades(:,1); % should be set to the column containing unix timestamps

%% Set variables
T = numel(times);
% Chop of times that don't match with size of trade window
chopOffTimes = rem(T,numTrades);
T = T-chopOffTimes;
times = times(chopOffTimes+1:end);
totIterations = T/numTrades; % must be an integer

% Zero out matrices
AIC = zeros(totIterations,1);
BIC = zeros(totIterations,1);
fullRMSE = zeros(totIterations,1);
mu =zeros(totIterations,1);
alpha = zeros(totIterations,1);
beta = zeros(totIterations,1);
expectedIntensity = zeros(totIterations,1);
branchingRatio = zeros(totIterations,1);
Nresponse = zeros(totIterations,1);
halfLife = zeros(totIterations,1);

%% Begin Iteration data points
beginTimes = 1;
endTimes = numTrades;
for iterationTrack = 1:totIterations
%% Initialize current loop varaibles
set = 1:numTrades;
timesNow = times(beginTimes:endTimes);
beginTimes = beginTimes + numTrades;
endTimes = endTimes + numTrades;

%% Introduce noise to matching time stamps
[~,uniqueStamps] = unique(timesNow);
idxNonunique = ~ismember(set,uniqueStamps);
for i = 1:numTrades
   if idxNonunique(i)
       timesNow(i) = timesNow(i) + rand;
   end
end
timesNow = sort(timesNow);

%% Call the minimization function
% set initial guess and reference function
mu0 = .5; % making mu >= or beta seems to work the best
alpha0 = .1; % keeping this about 1/10 of beta seems the best
beta0 = 1;
parameters = [mu0; alpha0; beta0];
func = @(parameters) HawkesMLE(parameters,timesNow);

%[fitParameters,logLikelihood,EXITFLAG] = fmincon(func,parameters,[-1 0 0],0,[],[],[],[],[],optimset('MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-8));
[fitParameters,logLikelihood,EXITFLAG] = fminunc(func,parameters,optimset('MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-8));

mu(iterationTrack) = fitParameters(1);
alpha(iterationTrack) = fitParameters(2);
beta(iterationTrack) = fitParameters(3);

% fmincon/fminunc appears to work the best

% Calculate the fitted conditional intensities
tempIntensity = zeros(numTrades,1);
for t = 2:numTrades
    tempIntensity(t) = exp(-beta(iterationTrack) * (timesNow(t) - timesNow(t-1)))*(1+tempIntensity(t-1));
end

% Calculate results
conditionalIntensity = mu(iterationTrack) + alpha(iterationTrack)*tempIntensity;
expectedIntensity(iterationTrack,1) = mu(iterationTrack)/(1-alpha(iterationTrack)/beta(iterationTrack));
branchingRatio(iterationTrack,1) = alpha(iterationTrack)/beta(iterationTrack);
Nresponse(iterationTrack,1) = alpha(iterationTrack)/(beta(iterationTrack)-alpha(iterationTrack));
halfLife(iterationTrack,1) = log(2)./ beta(iterationTrack);

%% Bin the data
% Add one for each empirical count
% Sum the conditional intensities within each window
empiricalBins = [];
fitBins = [];
binTimes = [];
binTimes(1) = timesNow(1);
timeTrack = timesNow(1);
empiricalBins(1) = 0;
fitBins(1) = 0;
n = 1;
for t = 1:numTrades
   if timesNow(t) < timeTrack
       empiricalBins(n) = empiricalBins(n) + 1;
       fitBins(n) = fitBins(n) + conditionalIntensity(t);
   else
       while timesNow(t) > timeTrack
       n = n + 1;
       timeTrack = timeTrack + 60;
       empiricalBins(n) = 0;
       binTimes(n) = timeTrack;
       fitBins(n) = 0;
       end
       empiricalBins(n) = 1;
       fitBins(n) = fitBins(n) + conditionalIntensity(t);
   end
end
numBins = 1:numel(empiricalBins);

%% Results
timeBegin{iterationTrack,1} = datestr(datenum([1970 1 1 0 0 timesNow(1)]));
timeEnd{iterationTrack,1} = datestr(datenum([1970 1 1 0 0 timesNow(end)]));
dateVector = [];
dateVector = zeros(numel(binTimes),1);
for i = 1:numel(binTimes)
    dateVector(i) = datenum([1970 1 1 0 0 binTimes(i)]);
end

setTitle = strcat('Begin:', timeBegin{iterationTrack},' End:', timeEnd{iterationTrack});

img(iterationTrack) = figure(iterationTrack);
%plot(dateVector,empiricalBins,dateVector,fitBins); % using actual dates caused distorition of the axes
plot(1:numel(dateVector),empiricalBins,1:numel(dateVector),fitBins)
title(setTitle)
xlabel('Time')
ylabel('Intensity')
legend('Empirical','Fit')
%datetick

[thisRMSE,thisAIC,thisBIC] = hawkesFitStatistics(empiricalBins,fitBins,logLikelihood);


fullRMSE(iterationTrack,1) = thisRMSE;
AIC(iterationTrack,1) = thisAIC;
BIC(iterationTrack,1) = thisBIC;

%% Save plot to image
if saveImages
    [~, finalMonth] = month(timeEnd);
    finalYear = year(timeEnd(iterationTrack));
    finalYear = num2str(finalYear);
    iter = num2str(iterationTrack);
    imageName = strcat('Image',iter,'_',finalMonth(iterationTrack,:),finalYear,'.jpg');
    figfile = fullfile(pathname,imageName);
    saveas(img(iterationTrack),figfile);
end
close all;
end

%% Output results
T = table(timeBegin,timeEnd,fullRMSE,AIC,BIC,mu,alpha,beta,expectedIntensity,branchingRatio,Nresponse,halfLife);
writetable(T,'Hawkes Output.csv');
figure
plot(fullRMSE);
title('RMSE Plot')
xlabel('Data Point')
ylabel('RMSE')

algoTime = toc %algoTime is the time in seconds for the entire script