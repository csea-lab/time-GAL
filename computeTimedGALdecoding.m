function [TimedGALdecoding] = computeTimedGALdecoding(condition1, condition2, trialsCondition1, trialsCondition2, channels, filename)
%   Compute Generalization Across Location
%   Generalization Across Location (GAL) describe the interrelation between
%   different locations of the brain showing similar patterns of
%   discrimination between two conditions.
%   
%   Input:
%       - subjectTest : index of the subject used as test data.
%       - condition1 : data matrix of channels*time*trials.
%       - condition2 : data matrix of channels*time*trials.
%       - trialsCondition1 : vector indicating who subject each trial
%       belongs to.
%       - trialsCondition1 : vector indicating who subject each trial
%       belongs to.
%       - channels : vector with indices of channels employed in analysis.
%       EXAMPLE: channels = [1:124 129];% Here, we use only 125 channels from Electrical Geodesics (EGI) 128-channel EEG
%
%   Number of trials does not need to match for each condition or subject.
%   Condtions A and B must contain the same number of time points and resolution.
%
%   Output:
%       - TimedGALdecoding
%           .GeneralizationMatrix : Generalization connectivity matrices
%           .CorrelationMatrix : Weight correlation matrix
%           .Topography : Decoding rate of each channel. Static and across
%           time
%


%%  1. Initialize 
    fprintf('\n 1. Initializing function. \n');

    % Start activating parallel pool for fast computation of subjects.
    parallelPool = gcp('nocreate');

    if not(isempty(parallelPool))
        fprintf('\n Parallel pool activated. Using %d workers ', parallelPool.NumWorkers);
    else
        fprintf('\n Parallel pool not activated.')
    end

    % Take register of time
    tic

%%  2. Compute GAL Generalization connectvity matrices

    fprintf('\n \n 2. Computing GAL matrices. Progress:\n');
    
    % Obtain number of total subjects,that should be contained on last
    % trial's subject.
    n = trialsCondition1(end); 
    
    % Pre-allocate output matrix
    generalizationConnectivityMatrices = zeros([n, channels(end), channels(end)]);

    parfor subjectTest = 1:n
        fprintf('\n Subject: %d out of %d   ', subjectTest, n);

        % Pre-allocate temporal matrix
        generalizationConnectivityMatrix = zeros([channels(end), channels(end)]); 
        
        % Obtain indices for training and test subject's data
        indTest1 = find(trialsCondition1 == subjectTest);
        indTest2 = find(trialsCondition2 == subjectTest);
        indTrain1 = setxor(1:length(trialsCondition1), indTest1);
        indTrain2 = setxor(1:length(trialsCondition2), indTest2);
    
        % Loop over all channels to train each one and test with the rest.
        for ch = channels 
           fprintf('.');
           % Train matrix corresponds with data from all subjects less the
           % selected subject.
            trainmatfull = [ transpose(squeeze(condition1(ch, :, indTrain1)));   transpose(squeeze(condition2(ch, :, indTrain2))) ];
           
           % Create a vector with labels 1 and 2 indicating belonging to
           % condition 1 or 2.
            labelvectorfull = [ ones([size(indTrain1), 1]) .* 1 ; ones([size(indTrain2), 1]) .* 2 ]; 
    
           % Train classifier object/model using Linear Discriminant Analysis (LDA). 
            obj = fitcdiscr(trainmatfull,labelvectorfull, "Prior", "empirical");
            
            fprintf('|');
            % Loop over all the channels to test the model and fill GAL.
            for tst = channels
                % Test matrix corresponds with data from the tested subject.
                testmatfull = [ transpose(squeeze(condition1(tst, :, indTest1)));   transpose(squeeze(condition2(tst, :, indTest2))) ];
                % Predict test information using trained object/model.
                Yhat = predict( obj, testmatfull );
                % Create a vector with true labels.
                Y = [ ones([size(indTest1), 1]) .* 1 ; ones([size(indTest2), 1]) .* 2 ];
                % Extract confusion matrix with true/predicted labels.
                aha = confusionmat(Y,Yhat);
                % Save the hit rate, calculated as correct conditions 1 or 2
                % divided by all cases.
                generalizationConnectivityMatrix(ch, tst) = sum(diag(aha))./(sum(sum(aha)));
            end
        end
        
        % Save result of each subject in same group of matrices
        generalizationConnectivityMatrices(subjectTest, :, :) = generalizationConnectivityMatrix;
    
    end


%%  3. Compute Correlation matrix of temporal contributions, i.e. weights
    
    fprintf('\n 3. Computing correlation matrix. Progress:\n');

    % Obtain the number of time points in data
    timebins = size(condition1,2)-1;
    
    % Pre-allocate space for temporal correlation results
    corr_vect = zeros([timebins, 1]);
    pcorr_vect = zeros([timebins, 1]);
    
    % Pre-allocate space for correlation results for all channels
    correlationR = zeros([length(channels), timebins]);
    correlationP = zeros([length(channels), timebins]);
    
    % Loop over all the channels to check their contribution(wiehgt) on time.
    for ch = channels
            fprintf('.');

            % Test matrix corresponds with data from all subjects.
            testmatfull = [ transpose(squeeze(condition1(ch, :, :)));   transpose(squeeze(condition2(ch, :, :))) ];
            
            % Vector indicating which condition is data.
            labelvectorfull2Corr = [  ones(size(condition1,3),1)*1;
                             ones(size(condition2,3),1)*2; ];

            % Loop over all time points to calculate contribution (weight)
            % of this time to classifier discrimination.
            for timebin = 1:timebins
                [r, p] = corrcoef(labelvectorfull2Corr, testmatfull(:,timebin));
                % Save temporal vectors of each time point.
                corr_vect(timebin) = r(2,1);
                pcorr_vect(timebin) = p(2,1);
            end

            % Save contribution (weight) of each electrode in each time
            % point along the trial.
            correlationR(ch, :) = corr_vect(:); % Pearson's correlation r
            correlationP(ch, :) = pcorr_vect(:); % Pearson's correlation p-value
    end


%%  4. Combining spatial and temporal information
    
    fprintf('\n 4. Combining spatial and temporal information. \n');

    % Pre-Allocate the channel*time decoding rate results
    topographyAcrossTime = zeros([channels(end), size(correlationR,2)]);

    % Loop over time to convolve the diagonal of decoding (hit rate of each
    % channel) in its temporal contribution.
    for t = 1:size(correlationR,2)
        topographyAcrossTime(:,t) = abs(correlationR(:,t)) .* diag(squeeze(mean(generalizationConnectivityMatrices)));
    end
    


%%  5. Finalize
    
    fprintf('\n 5. Finishing function. \n');

    % Deactivate parallel pool
    if not(isempty(parallelPool))
        delete(parallelPool)
    end
    

    % Wrap result in a same structure
    TimedGALdecoding = [];
    TimedGALdecoding.GeneralizationMatrix.GeneralizationConnectivityMatrices = generalizationConnectivityMatrices;
    TimedGALdecoding.CorrelationMatrix.CorrelationR = correlationR;
    TimedGALdecoding.CorrelationMatrix.CorrelationP = correlationP;
    TimedGALdecoding.Topography.TopographyStatic = size(diag(squeeze(mean(TimedGALdecoding.GeneralizationMatrix.GeneralizationConnectivityMatrices))));
    TimedGALdecoding.Topography.TopographyAcrossTime = topographyAcrossTime;

    % Show time elapsed
    elapsedtime = toc;
    fprintf('\n Computation of GAL matrices beloging to %d subjects finished. Time elapsed: %d minutes', n, elapsedtime/60)

    if nargin == 6
        fprintf(['\n Saving results in: ', filename, '\n'])
        save(filename, "TimedGALdecoding")
    end
end