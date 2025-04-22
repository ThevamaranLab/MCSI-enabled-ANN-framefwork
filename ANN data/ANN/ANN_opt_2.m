% ANN optimization (Shuffle-split procedure)

clc;
% Load input data (descriptors and mechanical property)
predictorData = Logdes; % Descriptors
targetData = LogSEA;    % Target property 

% Define the number of shuffles for Shuffle-Split validation
nShuffles = 1000;  % Number of shuffles

% Define the train-test split ratio
trainRatio = 0.8; % Approximately 64 points

% Define hyperparameters tuning range
numHiddenLayersList = 1:3;       % Number of hidden layers
numNeuronsList = 5:25;          % Number of neurons per layer
activationFunctions = {'tansig', 'logsig', 'poslin', 'purelin'}; % Different activation functions

% Store the performance metrics
results = struct();
results.Configurations = {};
results.TrainMSE = [];
results.TestMSE = [];
results.TrainRValue = [];
results.TestRValue = [];
results.TrainR2 = [];
results.TestR2 = [];
results.TrainMAPE = [];
results.TestMAPE = [];

% Grid search through hyper paramters combination
for numLayers = numHiddenLayersList
    for numNeurons = numNeuronsList
        for actFunc = activationFunctions
            trainMSEValues = zeros(1, nShuffles); % Store training MSE for each shuffle
            testMSEValues = zeros(1, nShuffles);  % Store testing MSE for each shuffle
            trainRValues = zeros(1, nShuffles);  % Store training R-value for each shuffle
            testRValues = zeros(1, nShuffles);   % Store testing R-value for each shuffle
            trainR2Values = zeros(1, nShuffles); % Store training R^2 for each shuffle
            testR2Values = zeros(1, nShuffles);  % Store testing R^2 for each shuffle
            trainMAPEValues = zeros(1, nShuffles); % Store training MAPE for each shuffle
            testMAPEValues = zeros(1, nShuffles);  % Store testing MAPE for each shuffle
            
            % Shuffle-Split cross-validation
            for shuffleIdx = 1:nShuffles
                % Randomly split the data
                indices = randperm(size(predictorData, 1));
                numTrain = round(trainRatio * size(predictorData, 1));
                trainIdx = indices(1:numTrain);
                testIdx = indices(numTrain+1:end);
                
                % Split data into training and testing sets
                trainX = predictorData(trainIdx, :)';
                trainY = targetData(trainIdx)';
                testX = predictorData(testIdx, :)';
                testY = targetData(testIdx)';
                
                % Construct the architecture
                hiddenLayerSizes = repmat(numNeurons, 1, numLayers);
                net = fitnet(hiddenLayerSizes, 'trainbr'); % Apply Bayesian regularization
                
                % Set activation functions for hidden layers
                for i = 1:numLayers
                    net.layers{i}.transferFcn = actFunc{1};
                end
                
                % Set training parameters
                net.divideParam.trainRatio = 0.8;
                net.divideParam.valRatio = 0.2;
                net.divideParam.testRatio = 0;
                net.trainParam.showWindow = false; % Suppress GUI
                
                % Train the network
                [net, ~] = train(net, trainX, trainY);
                
                % Evaluate the network on training set
                trainPredictions = net(trainX);
                trainMSEValues(shuffleIdx) = mean((trainPredictions - trainY).^2); % Training MSE
                trainRValues(shuffleIdx) = corr(trainPredictions', trainY');      % Training R-value
                trainR2Values(shuffleIdx) = 1 - sum((trainY - trainPredictions).^2) / sum((trainY - mean(trainY)).^2); % Training R^2
                trainMAPEValues(shuffleIdx) = mean(abs((trainY - trainPredictions) ./ trainY)) * 100; % Training MAPE
                
                % Evaluate the network on test set
                testPredictions = net(testX);
                testMSEValues(shuffleIdx) = mean((testPredictions - testY).^2);   % Testing MSE
                testRValues(shuffleIdx) = corr(testPredictions', testY');        % Testing R-value
                testR2Values(shuffleIdx) = 1 - sum((testY - testPredictions).^2) / sum((testY - mean(testY)).^2); % Testing R^2
                testMAPEValues(shuffleIdx) = mean(abs((testY - testPredictions) ./ testY)) * 100; % Testing MAPE
            end
            
            % Average metrics across all shuffles
            avgTrainMSE = mean(trainMSEValues);
            avgTestMSE = mean(testMSEValues);
            avgTrainRValue = mean(trainRValues);
            avgTestRValue = mean(testRValues);
            avgTrainR2 = mean(trainR2Values);
            avgTestR2 = mean(testR2Values);
            avgTrainMAPE = mean(trainMAPEValues);
            avgTestMAPE = mean(testMAPEValues);
            
            % Save the results
            results.Configurations{end+1} = sprintf('Layers: %d, Neurons: %d, ActFunc: %s', ...
                numLayers, numNeurons, actFunc{1});
            results.TrainMSE(end+1) = avgTrainMSE;
            results.TestMSE(end+1) = avgTestMSE;
            results.TrainRValue(end+1) = avgTrainRValue;
            results.TestRValue(end+1) = avgTestRValue;
            results.TrainR2(end+1) = avgTrainR2;
            results.TestR2(end+1) = avgTestR2;
            results.TrainMAPE(end+1) = avgTrainMAPE;
            results.TestMAPE(end+1) = avgTestMAPE;
            
            % Print metrics for the current configuration
            fprintf('Config: %s, Avg Train MSE: %.4f, Avg Test MSE: %.4f, Avg Train R: %.4f, Avg Test R: %.4f, Avg Train R^2: %.4f, Avg Test R^2: %.4f, Avg Train MAPE: %.2f%%, Avg Test MAPE: %.2f%%\n', ...
                results.Configurations{end}, avgTrainMSE, avgTestMSE, avgTrainRValue, avgTestRValue, avgTrainR2, avgTestR2, avgTrainMAPE, avgTestMAPE);
        end
    end
end

% Show the best configuration based on Test MSE
[~, bestIdx] = min(results.TestMSE);
fprintf('\nBest Configuration Based on Test MSE:\n');
fprintf('%s\nTrain MSE: %.4f, Test MSE: %.4f, Train R: %.4f, Test R: %.4f, Train R^2: %.4f, Test R^2: %.4f, Train MAPE: %.2f%%, Test MAPE: %.2f%%\n', ...
    results.Configurations{bestIdx}, results.TrainMSE(bestIdx), results.TestMSE(bestIdx), ...
    results.TrainRValue(bestIdx), results.TestRValue(bestIdx), ...
    results.TrainR2(bestIdx), results.TestR2(bestIdx), ...
    results.TrainMAPE(bestIdx), results.TestMAPE(bestIdx));
