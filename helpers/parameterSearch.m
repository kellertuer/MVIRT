function results = parameterSearch(startRanges, algorithmHandle, optValueHandle,inputValue,varargin)
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
% Manifold-Valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-12-21
% see LICENSE.txt

ip = inputParser;
addParameter(ip,'includeAllResults',false);
addParameter(ip,'MaxIterations',5);
addParameter(ip,'MinFunctionalValue',Inf);
addParameter(ip,'Export',[]);
parse(ip, varargin{:});
vars = ip.Results;
iterations = 0;
paramBoundaries = startRanges;
% Run for minimal values
resultsLeft.minParameters = paramBoundaries(1,:);
    debug('text',3,'Text',['--- Evaluating (lower bound) Parameters ',num2str(resultsLeft.minParameters),'...']);
resultsLeft.minOutput = algorithmHandle(resultsLeft.minParameters,inputValue);
resultsLeft.minValue = optValueHandle(resultsLeft.minParameters,resultsLeft.minOutput);
    debug('text',3,'Text',['--- ...yields ',num2str(resultsLeft.minValue)]);
% Run for maximal values.
resultsRight.minParameters = paramBoundaries(2,:);
    debug('text',3,'Text',['--- Evaluating (upper bound) Parameters ',num2str(resultsRight.minParameters),'...']);
resultsRight.minOutput = algorithmHandle(resultsRight.minParameters,inputValue);
resultsRight.minValue = optValueHandle(paramBoundaries(2,:),resultsRight.minOutput);
    debug('text',3,'Text',['--- ...yields ',num2str(resultsRight.minValue)]);

if resultsLeft.minValue < resultsRight.minValue
    results.minParameters = resultsLeft.minParameters;
    results.minValue = resultsLeft.minValue;
    results.minOutput = resultsLeft.minOutput;
else
    results.minParameters = resultsRight.minParameters;
    results.minValue = resultsRight.minValue;
    results.minOutput = resultsRight.minOutput;
end
if vars.includeAllResults
    results.allResults(1).Parameters = resultsRight.minParameters;
    results.allResults(1).Value = resultsRight.minValue;
    results.allResults(1).Output = resultsRight.minOutput;
    results.allResults(2).Parameters = resultsLeft.minParameters;
    results.allResults(2).Value = resultsLeft.minValue;
    results.allResults(2).Output = resultsLeft.minOutput;
end  
while (iterations < vars.MaxIterations) || (results.minValue>vars.MinFunctionalValue)
    iterations = iterations + 1;
    % evaluate mid point
    midParams = (resultsLeft.minParameters + resultsRight.minParameters)/2;
    debug('text',3,'Text',['--- Evaluating Parameters ',num2str(midParams),'...']);
    midOutput = algorithmHandle(midParams,inputValue);
    midValue = optValueHandle(midParams,midOutput);
    debug('text',3,'Text',['--- ...yields ',num2str(midValue)]);
    if resultsLeft.minValue < resultsRight.minValue % Left min -> omit right
        resultsRight.minParameters = midParams;
        resultsRight.minOutput = midOutput;
        resultsRight.minValue = midValue;
        % save min
        results.minParameters = resultsLeft.minParameters;
        results.minOutput = resultsLeft.minOutput;
        results.minValue = resultsLeft.minValue;
    else %if resultsLeft.minValue > resultsRight.minValue % right min -> omit left
        resultsLeft.minParameters = midParams;
        resultsLeft.minOutput = midOutput;
        resultsLeft.minValue = midValue;
        % save min
        results.minParameters = resultsRight.minParameters;
        results.minOutput = resultsRight.minOutput;
        results.minValue = resultsRight.minValue;
    end
    if vars.includeAllResults
        results.allResults(iterations+2).Parameters = midParams;
        results.allResults(iterations+2).Value = midValue;
        results.allResults(iterations+2).Output = midOutput;
    end
    if ~isempty(vars.Export)
        vars.Export(['-',num2str(iterations)],midOutput);
    end
end
    if midValue < results.minValue % midpoint from last iteration even smaller -> save here, 
        results.minParameters = midParams;
        results.minOutput = midOutput;
        results.minValue = midValue;
    end
    debug('text',3,'Text',['The resulting minimum (',num2str(results.minValue),') is obtained for the parameters ',num2str(results.minParameters),'.']);
    if ~isempty(vars.Export)
        vars.Export('-minimum',results.minOutput);
    end
end
