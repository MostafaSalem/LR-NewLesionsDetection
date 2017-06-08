function [trainedClassifier, validationAccuracy] = trainClassifier3(trainingData)

inputTable = trainingData;
predictorNames = {'T1_Basal', 'T2_Basal', 'PD_Basal', 'FLAIR_Basal',...
                  'T1_FollowUp', 'T2_FollowUp', 'PD_FollowUp', 'FLAIR_FollowUp',...
                  'T1_Diff', 'T2_Diff', 'PD_Diff', 'FLAIR_Diff',...
                  'T1_Jacobian','T2_Jacobian','PD_Jacobian','FLAIR_Jacobian',...
                  'T1_Divergence','T2_Divergence','PD_Divergence','FLAIR_Divergence',...
                  'T1_NormDiv','T2_NormDiv','PD_NormDiv','FLAIR_NormDiv'};

predictors = inputTable(:, predictorNames);
response = inputTable.Response;
isCategoricalPredictor = [false, false, false, false, false, false, false, false,false,false,false, false, false, false, false, false, false, false,false,false,false, false, false, false];


successClass = '1';
failureClass = '0';
missingClass = '';
responseCategories = categories(response);
successFailureAndMissingClasses = categorical({failureClass; successClass; missingClass}, responseCategories);
isMissing = isundefined(response);
zeroOneResponse = double(ismember(response, successClass));
zeroOneResponse(isMissing) = NaN;
% Prepare input arguments to fitglm.
categoricalPredictorIndex = find(isCategoricalPredictor);
concatenatedPredictorsAndResponse = [predictors, table(zeroOneResponse)];

% Train using zero-one responses, specifying which predictors are
% categorical.
GeneralizedLinearModel = fitglm(...
    concatenatedPredictorsAndResponse, ...
    'Distribution', 'binomial', ...
    'link', 'logit', ...
    'CategoricalVars', categoricalPredictorIndex);

% Convert predicted probabilities to predicted class labels and scores.
convertSuccessProbsToPredictions = @(p) successFailureAndMissingClasses( ~isnan(p).*( (p>=0.5) + 1 ) + isnan(p)*3 );
returnMultipleValuesFcn = @(varargin) varargin{1:max(1,nargout)};
scoresFcn = @(p) [p, 1-p];
predictionsAndScoresFcn = @(p) returnMultipleValuesFcn( convertSuccessProbsToPredictions(p), scoresFcn(p) );

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
logisticRegressionPredictFcn = @(x) predictionsAndScoresFcn( predict(GeneralizedLinearModel, x) );
trainedClassifier.predictFcn = @(x) logisticRegressionPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'T1_Basal', 'T2_Basal', 'PD_Basal', 'FLAIR_Basal',...
                  'T1_FollowUp', 'T2_FollowUp', 'PD_FollowUp', 'FLAIR_FollowUp',...
                  'T1_Diff', 'T2_Diff', 'PD_Diff', 'FLAIR_Diff',...
                  'T1_Jacobian','T2_Jacobian','PD_Jacobian','FLAIR_Jacobian',...
                  'T1_Divergence','T2_Divergence','PD_Divergence','FLAIR_Divergence',...
                  'T1_NormDiv','T2_NormDiv','PD_NormDiv','FLAIR_NormDiv'};
trainedClassifier.GeneralizedLinearModel = GeneralizedLinearModel;
trainedClassifier.SuccessClass = successClass;
trainedClassifier.FailureClass = failureClass;
trainedClassifier.MissingClass = missingClass;
trainedClassifier.ClassNames = {successClass; failureClass};
trainedClassifier.About = 'This struct is a trained classifier exported from Classification Learner R2016a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedClassifier''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');
validationAccuracy=1;

% Extract predictors and response
% This code processes the data into the right shape for training the
% classifier.
% inputTable = trainingData;
% predictorNames = {'T1_Basal','T2_Basal','PD_Basal','FLAIR_Basal','T1_FollowUp', 'T2_FollowUp', 'PD_FollowUp', 'FLAIR_FollowUp', 'T1_Diff', 'T2_Diff', 'PD_Diff', 'FLAIR_Diff'};
% predictors = inputTable(:, predictorNames);
% response = inputTable.Response;
% isCategoricalPredictor = [false, false, false, false,false, false, false, false, false, false, false, false];

% % Perform cross-validation
% KFolds = 2;%size(inputTable,1);
% cvp = cvpartition(response, 'KFold', KFolds);
% % Initialize the predictions and scores to the proper sizes
% validationPredictions = response;
% numObservations = size(predictors, 1);
% numClasses = 2;
% validationScores = NaN(numObservations, numClasses);
% for fold = 1:KFolds
%     trainingPredictors = predictors(cvp.training(fold), :);
%     trainingResponse = response(cvp.training(fold), :);
%     foldIsCategoricalPredictor = isCategoricalPredictor;
%     
%     % Train a classifier
%     % This code specifies all the classifier options and trains the classifier.
%     % For logistic regression, the response values must be converted to zeros
%     % and ones because the responses are assumed to follow a binomial
%     % distribution.
%     % 1 or true = 'successful' class
%     % 0 or false = 'failure' class
%     % NaN - missing response.
%     successClass = '1';
%     failureClass = '0';
%     missingClass = '';
%     responseCategories = categories(trainingResponse);
%     successFailureAndMissingClasses = categorical({failureClass; successClass; missingClass}, responseCategories);
%     isMissing = isundefined(trainingResponse);
%     zeroOneResponse = double(ismember(trainingResponse, successClass));
%     zeroOneResponse(isMissing) = NaN;
%     % Prepare input arguments to fitglm.
%     categoricalPredictorIndex = find(foldIsCategoricalPredictor);
%     concatenatedPredictorsAndResponse = [trainingPredictors, table(zeroOneResponse)];
%     
%     % Train using zero-one responses, specifying which predictors are
%     % categorical.
%     GeneralizedLinearModel = fitglm(...
%         concatenatedPredictorsAndResponse, ...
%         'Distribution', 'binomial', ...
%         'link', 'logit', ...
%         'CategoricalVars', categoricalPredictorIndex);
%     
%     % Convert predicted probabilities to predicted class labels and scores.
%     convertSuccessProbsToPredictions = @(p) successFailureAndMissingClasses( ~isnan(p).*( (p>=0.5) + 1 ) + isnan(p)*3 );
%     returnMultipleValuesFcn = @(varargin) varargin{1:max(1,nargout)};
%     scoresFcn = @(p) [p, 1-p];
%     predictionsAndScoresFcn = @(p) returnMultipleValuesFcn( convertSuccessProbsToPredictions(p), scoresFcn(p) );
%     
%     % Create the result struct with predict function
%     logisticRegressionPredictFcn = @(x) predictionsAndScoresFcn( predict(GeneralizedLinearModel, x) );
%     validationPredictFcn = @(x) logisticRegressionPredictFcn(x);
%     
%     % Add additional fields to the result struct
%     
%     % Compute validation predictions and scores
%     validationPredictors = predictors(cvp.test(fold), :);
%     [foldPredictions, foldScores] = validationPredictFcn(validationPredictors);
%     
%     % Store predictions and scores in the original order
%     validationPredictions(cvp.test(fold), :) = foldPredictions;
%     validationScores(cvp.test(fold), :) = foldScores;
% end
% 
% correctPredictions = (validationPredictions == response);
% validationAccuracy = sum(correctPredictions)/length(correctPredictions);
