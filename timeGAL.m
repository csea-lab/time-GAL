function [timeGALoutput] = timeGAL(condition1, condition2, trialsCondition1, trialsCondition2, varargin)
%   Compute Generalization Across Location
%   Generalization Across Location (GAL) describe the interrelation between
%   different locations of the brain showing similar patterns of
%   discrimination between two conditions.
%   
%   Input:
%       - condition1 : data matrix of channels*time*trials.
%       - condition2 : data matrix of channels*time*trials.
%       - trialsCondition1 : vector indicating who subject each trial
%       belongs to.
%       - trialsCondition1 : vector indicating who subject each trial
%       belongs to.
%
%   Parameters:
%           - 'ParallelComputing' = true/false.
%           - 'ParallelComputingCores' = number of workers/cores to create parpool.
%           - 'Channels' = vector that indicates the computation of only certain channels
%                   Example: channels = [1:124 129];% Here, we use only 125 channels from Electrical Geodesics (EGI) 128-channel EEG
%           - 'Time' = vector that indicates the temporal windows of analysis.
%                   Example: time = [300:500];% Here, we use use 200 time points
%           - 'FileName' = file name where save output of the TimeGAL toolbox.
%           - 'AlphaGAL' = threshold for statistical analysis of GAL matrices (default = 0.05), corrected by Bonferroni.
%           - 'AlphaCorrelation' = threshold for statistical analysis of Correlation matrices (default = 0.01)
%
%   Output:
%        GeneralizationMatrix.GAL = GAL matrices (subject by channel by channel).
%        GeneralizationMatrix.Topography = Diagonal of GAL.
%        GeneralizationMatrix.GALmask = Mask of statistical significant GAL.
%        GeneralizationMatrix.GALmaskPos =  Mask of direct statisticallly significant GAL.
%        GeneralizationMatrix.GALmaskNeg = Mask of inverse statistically significant GAL.
%        
%        CorrelationMatrix.CorrelationR = Pearson's correlation (channel by time).
%        CorrelationMatrix.CorrelationP = P-value of correlation (channel by time).
%        CorrelationMatrix.CorrelationMask = Mask of statistically significant correlations.
%        
%        TimeGAL.TimeGAL = Time-GAL decoding (channel by channel by time).
%        TimeGAL.TimeTopography = Topography across time (channel by time).
%        
%        Parameters.alphaGAL
%        Parameters.alphaCorrelation
%        Parameters.Parallel 
%        Parameters.ParallelCores
%        Parameters.Channels
%        Parameters.Time
%        Parameters.FileName
%        Parameters.ListOfSubjects = List of which subjects where included.
%        Parameters.InputAverage = [{mean(condition1,3)},{mean(condition2,3)}],
%
%
%   Number of trials does not need to match for each condition or subject.
%   However, huge difference between amount of trials between conditions
%   can lead to classifier bias, resulting in ficticious high decoding
%   rates.
%   trialsCondition list of subjects can use any number for each subject
%   and it is not restricted to start from 1 nor increase by one unit.
%   Condtions A and B must contain the same number of time points and resolution.
%
%   Example:
%   
%   [timeGALoutput] = timeGAL(condition1, condition2, trialsCondition1, trialsCondition2, 'ParallelComputing', true, 'ParallelComputingCores', 4, 'Channels', [1:124 129], 'Filename', 'resultsTimeGAL.mat')



%% Arguments

    arg = inputParser;
    addParameter(arg, 'Channels', [], @(x) isnumeric(x) && isvector(x)); % 
    addParameter(arg, 'Time', [], @(x) isnumeric(x) && isvector(x)); % 
    addParameter(arg, 'FileName', '', @ischar); % 
    addParameter(arg, 'AlphaGAL', 0.05, @numeric); % 
    addParameter(arg, 'AlphaCorrelation', 0.01, @numeric); % 
    addParameter(arg, 'ParallelComputing', false, @islogical); % 
    addParameter(arg, 'ParallelComputingCores', [], @isnumeric); % 


    % Access to arguments
    parse(arg, varargin{:});
    parallel = arg.Results.ParallelComputing;
    parallelCores = arg.Results.ParallelComputingCores;
    if isempty(arg.Results.Channels)
        channels = 1:size(condition1, 1); 
    else
        channels = arg.Results.Channels; 
    end
    if ~isempty(arg.Results.Time)
        condition1 = condition1(:,arg.Results.Time,:); 
        condition2 = condition2(:,arg.Results.Time,:);
        timeWindow = arg.Results.Time;
    else
        timeWindow = 1:size(condition1,2);
    end
    filename = arg.Results.FileName;
    alphaGAL = arg.Results.AlphaGAL;
    alphaCorrelation = arg.Results.AlphaCorrelation;

    if length(trialsCondition1) > length(trialsCondition2)*2 || length(trialsCondition2) > length(trialsCondition1)*2
        warning(['TimeGAL alert: One of the conditions containes more than double of the observations contained in the other condition.' ...
            ' Big differences in amount of trials between conditions can lead to classifier bias, resulting in ficticious high decoding rates.'])
    end
    


%%  1. Initialize 
    fprintf('\n 1. Initializing function. \n');

    % Take register of time
    t1 = tic;

    % Verificar si el Parallel Computing Toolbox está disponible
    isParallelToolboxAvailable = license('test', 'Distrib_Computing_Toolbox');

    if parallel && isParallelToolboxAvailable
        fprintf('\n Using Parallel Computing Toolbox with %d cores \n\n', parallelCores);
        pool = gcp('nocreate');
        if isempty(pool)
            parpool(parallelCores);
        elseif pool.NumWorkers ~= parallelCores
            delete(pool);
            parpool(parallelCores);
        end
    else
        pool = gcp('nocreate');
        if ~isempty(pool)
            delete(pool);
        end
    end


%%  2. Compute GAL Generalization connectvity matrices

    fprintf('\n \n 2. Computing GAL matrices. Progress:\n');
    
    % Obtain number of total subjects,that should be contained on last
    % trial's subject.
    listOfSubjects = unique(trialsCondition1);
    n = length(listOfSubjects); 
    
    % Pre-allocate output matrix
    generalizationConnectivityMatrices = zeros([n, channels(end), channels(end)]);
    
    if parallel % PARALLEL 
        parfor subjectTest = 1:n
            fprintf('\n Subject: %d | %d out of %d  \n', listOfSubjects(subjectTest), subjectTest, n);
            t2 = tic;
    
            % Pre-allocate temporal matrix
            generalizationConnectivityMatrix = zeros([channels(end), channels(end)]); 
            
            % Obtain indices for training and test subject's data
            indTest1 = find(trialsCondition1 == listOfSubjects(subjectTest));
            indTest2 = find(trialsCondition2 == listOfSubjects(subjectTest));
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
                labelvectorfull = [ ones([length(indTrain1), 1]) .* 1 ; ones([length(indTrain2), 1]) .* 2 ]; 

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
            subTimeElapsed = toc(t2)
            fprintf('\n Time elapse: %d minutes %d seconds', floor(subTimeElapsed / 60), mod(subTimeElapsed, 60));
        end
    else % NO PARALLEL 
        for subjectTest = 1:n
            fprintf('\n Subject: %d | %d out of %d  \n', listOfSubjects(subjectTest), subjectTest, n);
            t2 = tic;
    
            % Pre-allocate temporal matrix
            generalizationConnectivityMatrix = zeros([channels(end), channels(end)]); 
            
            % Obtain indices for training and test subject's data
            indTest1 = find(trialsCondition1 == listOfSubjects(subjectTest));
            indTest2 = find(trialsCondition2 == listOfSubjects(subjectTest));
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
                labelvectorfull = [ ones([length(indTrain1), 1]) .* 1 ; ones([length(indTrain2), 1]) .* 2 ]; 
        
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
                    Y = [ ones([length(indTest1), 1]) .* 1 ; ones([length(indTest2), 1]) .* 2 ];
                    % Extract confusion matrix with true/predicted labels.
                    aha = confusionmat(Y,Yhat);
                    % Save the hit rate, calculated as correct conditions 1 or 2
                    % divided by all cases.
                    generalizationConnectivityMatrix(ch, tst) = sum(diag(aha))./(sum(sum(aha)));
                end
            end
            
            % Save result of each subject in same group of matrices
            generalizationConnectivityMatrices(subjectTest, :, :) = generalizationConnectivityMatrix;
            subTimeElapsed = toc(t2)
            fprintf('\n Time elapse: %d minutes %d seconds', floor(subTimeElapsed / 60), mod(subTimeElapsed, 60));
        end
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
    

%%  4. Combining spatial and temporal information
    
    fprintf('\n 5. Calculating statistics and output data. \n');
    
    
    chanceLevel = 0.5;
    
    % T-Test contrast against chance Level (0.5)
    alpha_corrected = alphaGAL/125; % Bonferroni correction.
    [~,p] = ttest(generalizationConnectivityMatrices(:,channels,channels), chanceLevel, 'tail', 'right');
    maskpos = squeeze(p < alpha_corrected);
    [~,p] = ttest(generalizationConnectivityMatrices(:,channels,channels), chanceLevel, 'tail', 'left');
    maskneg = (squeeze(p < alpha_corrected)  );
    maskposneg = maskpos + ( maskneg  * - 1 );

    % Calculation of Time-GAL
    maskGAL = abs(maskposneg);%(channels,channels));
    maskTime = correlationP < alphaCorrelation;
    %%indcs = find(sum(maskGAL)>=1);
    [row, col] = find(maskGAL);
    TimeGAL = zeros(125,125,size(maskTime,2));
    for i = 1:length(row)
        TimeGAL(row(i),col(i), maskTime(row(i),:) == 1) = 1;
    end
    TimeGAL = TimeGAL .* maskposneg; %(channels,channels);

    
    
%%  5. Finalize
    
    fprintf('\n 6. Finishing function. \n');

    % Deactivate parallel pool
    if parallel && isParallelToolboxAvailable
        if ~isempty(pool)
            delete(pool);
        end
    end
    

    % Wrap result in a same structure
    timeGALoutput = [];
    timeGALoutput.GeneralizationMatrix.GAL = generalizationConnectivityMatrices;
    timeGALoutput.GeneralizationMatrix.Topography = diag(squeeze(mean(generalizationConnectivityMatrices)));
    timeGALoutput.GeneralizationMatrix.GALmask = maskposneg;
    timeGALoutput.GeneralizationMatrix.GALmaskPos = maskpos;
    timeGALoutput.GeneralizationMatrix.GALmaskNeg = maskneg;

    timeGALoutput.CorrelationMatrix.CorrelationR = correlationR;
    timeGALoutput.CorrelationMatrix.CorrelationP = correlationP;
    timeGALoutput.CorrelationMatrix.CorrelationMask = maskTime;

    timeGALoutput.TimeGAL.TimeGAL = TimeGAL;
    timeGALoutput.TimeGAL.TimeTopography = topographyAcrossTime;

    timeGALoutput.Parameters.alphaGAL = alphaGAL;
    timeGALoutput.Parameters.alphaCorrelation = alphaCorrelation;
    timeGALoutput.Parameters.Parallel = parallel;
    timeGALoutput.Parameters.ParallelCores = parallelCores;
    timeGALoutput.Parameters.Channels = channels;
    timeGALoutput.Parameters.Time = timeWindow;
    timeGALoutput.Parameters.FileName = filename;
    timeGALoutput.Parameters.ListOfSubjects = listOfSubjects;
    timeGALoutput.Parameters.InputAverage = [{mean(condition1,3)},{mean(condition2,3)}];

    % Show time elapsed
    elapsedTime = toc(t1);
    fprintf('\n Computation of GAL matrices beloging to %d subjects finished. Time elapsed: %d minutes and % seconds', n, floor(elapsedTime / 60), mod(elapsedTime, 60))

    if ~isempty(filename)
        fprintf(['\n Saving results in: ', filename, '\n'])
        save(filename, "timeGALoutput")
    end
    
    % Final warnings
    if n < 20
        warning(['TimeGAL alert: The sample is smaller than 20 subjects. TimeGAL toolbox use the parametric T-Test contrast to calculate wether each decoding rate is higher than chance level.' ...
            ' However, user should ensure statistical assumptions for the T-Test or calculate non-parametrical alternatives using the GAL matrices.'])
    end
    if length(trialsCondition1) > length(trialsCondition2)*2 || length(trialsCondition2) > length(trialsCondition1)*2
        warning(['TimeGAL alert: One of the conditions containes more than double of the observations contained in the other condition.' ...
            ' Big differences in amount of trials between conditions can lead to classifier bias, resulting in ficticious high decoding rates.'])
    end

end