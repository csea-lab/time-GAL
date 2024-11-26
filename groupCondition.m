function [conditionData, listSubject] = groupCondition(path, condition, matname, varargin)
% Read files from path and condition to create one file 
%
% Input:
%   path = './'
%   condition = 'Data_Pleasant'
%   filename = 'ConditionPleasant'
% Output:
%   Condition [channel,time,trials]
%   condition = 'Data_Pleasant'
%
%  [conditionData, listSubject] = groupCondition(path, condition, matname, varargin)

arg = inputParser;
addParameter(arg, 'FileName', '', @ischar); %
addParameter(arg, 'ResampleProportion', [], @isnumeric); %
addParameter(arg, 'TimeWindow', [], @isnumeric); %

parse(arg, varargin{:});
filename = arg.Results.FileName;
resampleProportion = arg.Results.ResampleProportion;
timeWin = arg.Results.TimeWindow;




files = dir([path,'*',condition, '*.mat']);
fprintf(['\n \n A total of %d files were found with tag: ', condition], length(files))

tmp = load(strcat(files(1).folder, '/',files(1).name));
disp(tmp)
eval(strcat(['nchannels = size(tmp.', matname, ', 1);']))
eval(['ntime = size(tmp.', matname, ', 2);'])

conditionData = zeros([nchannels,ntime, 1]);
listSubject = zeros([1, 1]);

for f = 1:length(files)
    fprintf(['\n Reading file: ', files(f).name, ' | File %d out of %d'], f, length(files))
    tmp = load(strcat(files(f).folder, '/',files(f).name));
    eval(['tmpData = tmp.', matname, ';'])
    
    if ~isempty(timeWin) 
        tmpData = tmpData(:,timeWin,:);
            if ~isempty(resampleProportion)
                [p, q] = rat(resampleProportion);
                tmpData = resample(tmpData, p, q, 'Dimension',2);
            end
    end

    conditionData = cat(3, conditionData, tmpData);
    listSubject = cat(1, listSubject, ones(size(tmpData, 3), 1) .* f);
end
conditionData(:,:,1) = []; listSubject(1) = [];

if ~isempty(filename)
    eval(['Data_', filename,' = conditionData;'])
    eval(['List_', filename,' = listSubject;'])
    fprintf(['\n \n Saving conditionData matrix and listSubjects vector as file: ', filename])
    save(filename, ['Data_', filename], ['List_', filename])
end

