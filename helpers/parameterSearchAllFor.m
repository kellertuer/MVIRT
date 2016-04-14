function results = parameterSearchAllFor(startRanges, algorithmHandle, optValueHandle,inputValue,varargin)
% parameterCheck(startRanges, functional) - performs a parameter Check by
% binary division of the intervals for the specified parameters in
% startRanges and calles the algorithm and optimizes with respect to
% 
% INPUT
%   startRanges : Parameters that are to be checked simultaneously
%        in the form of columns having the minimal and maximal starting
%        value: [ a_min, b_min, c_min; a_max, b_max, c_max
%   algorithm   : the algorithm that is run - the function handle has to be
%       of the form resultAlg = algorithm(paramVec,inputValue)
%   functional  : functional, that takes the parameters and the result of
%       of algorithm to determine the optimality - the functional has to be
%        of the form nonnegativeValue = functional(paramVec,resultAlg)
%   inputValue  : data to perform the algorithm on
%
% OUTPUT
%   results: a struct with
%       results.minParameters containing the minimal parameters
%       results.minValue containing the corresponding functional value
%       results.minOutput containing the result of the algorithm for the
%           optimal parameters
%
% OPTIONAL PARAMETERS
%   'includeAllResults' : (false) include all results in the struct using
%       results.allResults - a structure erray consisting of the fields
%           Parameters, Value, and Output, i.e. each entry is structured
%           as the results struct (see above)
%   'Export': ([]) a function handle accepting 2 parameters, a short
%          string (iteration or result) and the result and produces a
%          certain type of export - a figure or a file or something like that
%   'MaxIterations': (5) number of recursive splittings of the parameter
%       ranges as a stopping criterion
%   'MinFunctionalValue': (Inf) lower bound for the functional value to be
%   reached as a stopping criterion
% If both stopping criteria are set, this function stops if one of the criteria
% is met.
%
% ---
% ManImRes 1.0
% R. Bergmann ~ 2015-12-21

ip = inputParser;
addParameter(ip,'includeAllResults',false);
addParameter(ip,'MaxIterations',5);
addParameter(ip,'MinFunctionalValue',Inf);
addParameter(ip,'Export',[]);
parse(ip, varargin{:});
vars = ip.Results;
paramBoundaries = startRanges;
% Run for minimal values
resultsLeft.minParameters = paramBoundaries(1,:);
    debug('text',1,'Text',['--- Evaluating (lower bound) Parameters ',num2str(resultsLeft.minParameters),'...']);
resultsLeft.minOutput = algorithmHandle(resultsLeft.minParameters,inputValue);
resultsLeft.minValue = optValueHandle(resultsLeft.minParameters,resultsLeft.minOutput);
    debug('text',1,'Text',['--- ...yields ',num2str(resultsLeft.minValue)]);
% Run for maximal values.
resultsRight.minParameters = paramBoundaries(2,:);
    debug('text',1,'Text',['--- Evaluating (upper bound) Parameters ',num2str(resultsRight.minParameters),'...']);
resultsRight.minOutput = algorithmHandle(resultsRight.minParameters,inputValue);
resultsRight.minValue = optValueHandle(paramBoundaries(2,:),resultsRight.minOutput);
    debug('text',1,'Text',['--- ...yields ',num2str(resultsRight.minValue)]);

if resultsLeft.minValue < resultsRight.minValue
    results.minParam = resultsLeft.minParameters;
    results.minValue = resultsLeft.minValue;
    results.minOutput = resultsLeft.minOutput;
else
    results.minParam = resultsRight.minParameters;
    results.minValue = resultsRight.minValue;
    results.minOutput = resultsRight.minOutput;
end
%if vars.includeAllResults
    results.allResults(1).Parameters = resultsRight.minParameters;
    results.allResults(1).Value = resultsRight.minValue;
    results.allResults(1).Output = resultsRight.minOutput;
    results.allResults(2).Parameters = resultsLeft.minParameters;
    results.allResults(2).Value = resultsLeft.minValue;
    results.allResults(2).Output = resultsLeft.minOutput;
%end
leftB = resultsLeft.minParameters;
rightB = resultsRight.minParameters;
pos=1;
iterations=0;
results.minValue = Inf;
while (iterations < vars.MaxIterations) || (results.minValue>vars.MinFunctionalValue)
    % run through all mid points look for the best
    paramChoices = repmat({[0 1/4 1/2 3/4 1]},1,length(leftB));
    c = cell(1, numel(paramChoices));
    [c{:}] = ndgrid( paramChoices{:} );
    paramChoices = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
    bestMidValue = Inf;
    bestMidOutput = [];
    iterations=iterations+1;
    sizeV = rightB-leftB;
    for i=1:size(paramChoices,1)
%        paramPermutations = perms(paramChoices(i,:));
%        for j=1:size(paramPermutations,1)
%           midParams = leftB + paramPermutations(j,:).*sizeV;
            midParams = leftB + paramChoices(i,:).*sizeV;
            found=false;
            for k=1:length(results.allResults)
                if all(results.allResults(k).Parameters== midParams);
                debug('text',1,'Text',['--- Parameters ',num2str(midParams),' already evaluated.']);
                    found=true;
                    midOutput = results.allResults(k).Output;
                    midValue = results.allResults(k).Value;
                    break;
                end
            end
            if ~found % compute
                debug('text',1,'Text',['--- Evaluating Parameters ',num2str(midParams),'...']);
                midOutput = algorithmHandle(midParams,inputValue);
                midValue = optValueHandle(midParams,midOutput);
                results.allResults(pos+2).Parameters = midParams;
                results.allResults(pos+2).Value = midValue;
                results.allResults(pos+2).Output = midOutput;
                pos = pos + 1;
            end
            if midValue < bestMidValue
                bestMidParams = midParams;
                bestMidOutput = midOutput;
                bestMidValue = midValue;
            end
            if ~isempty(vars.Export)
                vars.Export(['-',num2str(pos)],midOutput);
            end
            debug('text',1,'Text',['--- ...yields ',num2str(midValue)]);
%       end
    end
    leftB = max(bestMidParams-sizeV./4,leftB);
    rightB = min(leftB + sizeV/2, rightB);
    
    % Having the optimal inner value -> save, update boundaries
    if bestMidValue < results.minValue
        results.minParam = bestMidParams;
        results.minValue = bestMidValue;
        results.minOutput = bestMidOutput;
    end %otherwise: parameter range is smaller, lets start new
end
%if ~vars.includeAllResults
%    rmfield(results,'allResults');
%end
    debug('text',1,'Text',['The resulting minimum (',num2str(results.minValue),') is obtained for the parameters ',num2str(results.minParam),'.']);
    if ~isempty(vars.Export)
        vars.Export('-minimum',results.minOutput);
    end
end