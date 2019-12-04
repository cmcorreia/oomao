function b = padarray(varargin)
%PADARRAY Pad an array.
%   B = PADARRAY(A,PADSIZE) pads array A with PADSIZE(k) number of zeros
%   along the k-th dimension of A.  PADSIZE should be a vector of
%   positive integers.
%
%   B = PADARRAY(A,PADSIZE,PADVAL) pads array A with PADVAL (a scalar)
%   instead of with zeros.
%
%   B = PADARRAY(A,PADSIZE,PADVAL,DIRECTION) pads A in the direction
%   specified by the string DIRECTION.  DIRECTION can be one of the
%   following strings.  
%
%       String values for DIRECTION
%       'pre'         Pads before the first array element along each
%                     dimension .
%       'post'        Pads after the last array element along each
%                     dimension. 
%       'both'        Pads before the first array element and after the
%                     last array element along each dimension.
%
%   By default, DIRECTION is 'both'.
%
%   B = PADARRAY(A,PADSIZE,METHOD,DIRECTION) pads array A using the
%   specified METHOD.  METHOD can be one of these strings:
%
%       String values for METHOD
%       'circular'    Pads with circular repetion of elements.
%       'replicate'   Repeats border elements of A.
%       'symmetric'   Pads array with mirror reflections of itself. 
% 
%   Class Support
%   -------------
%   When padding with a constant value, A can be numeric or logical.
%   When padding using the 'circular', 'replicate', or 'symmetric'
%   methods, A can be of any class.  B is of the same class as A.
%
%   Example
%   -------
%   Add three elements of padding to the beginning of a vector.  The
%   padding elements contain mirror copies of the array.
%
%       b = padarray([1 2 3 4],3,'symmetric','pre')
%
%   Add three elements of padding to the end of the first dimension of
%   the array and two elements of padding to the end of the second
%   dimension.  Use the value of the last array element as the padding
%   value.
%
%       B = padarray([1 2; 3 4],[3 2],'replicate','post')
%
%   Add three elements of padding to each dimension of a
%   three-dimensional array.  Each pad element contains the value 0.
%
%       A = [1 2; 3 4];
%       B = [5 6; 7 8];
%       C = cat(3,A,B)
%       D = padarray(C,[3 3],0,'both')
%
%   See also CIRCSHIFT, IMFILTER.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.11.4.3 $  $Date: 2003/08/23 05:53:08 $

[a, method, padSize, padVal, direction] = ParseInputs(varargin{:});

if isempty(a),% treat empty matrix similar for any method

   if strcmp(direction,'both')
      sizeB = size(a) + 2*padSize;
   else
      sizeB = size(a) + padSize;
   end

   b = mkconstarray(class(a), padVal, sizeB);
   
else
  switch method
    case 'constant'
        b = ConstantPad(a, padSize, padVal, direction);
        
    case 'circular'
        b = CircularPad(a, padSize, direction);
  
    case 'symmetric'
        b = SymmetricPad(a, padSize, direction);
        
    case 'replicate'
        b = ReplicatePad(a, padSize, direction);
  end      
end

if (islogical(a))
    b = logical(b);
end

%%%
%%% ConstantPad
%%%
function b = ConstantPad(a, padSize, padVal, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
sizeB = zeros(1,numDims);
for k = 1:numDims
    M = size(a,k);
    switch direction
        case 'pre'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + padSize(k);
            
        case 'post'
            idx{k}   = 1:M;
            sizeB(k) = M + padSize(k);
            
        case 'both'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + 2*padSize(k);
    end
end

% Initialize output array with the padding value.  Make sure the
% output array is the same type as the input.
b         = mkconstarray(class(a), padVal, sizeB);
b(idx{:}) = a;


function out = mkconstarray(class, value, size)
%MKCONSTARRAY creates a constant array of a specified numeric class.
%   A = MKCONSTARRAY(CLASS, VALUE, SIZE) creates a constant array 
%   of value VALUE and of size SIZE.

%   Copyright 1993-2003 The MathWorks, Inc.  
%   $Revision: 1.8.4.1 $  $Date: 2003/01/26 06:00:35 $

out = repmat(feval(class, value), size);

%%%
%%% CircularPad
%%%
function b = CircularPad(a, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  dimNums = 1:M;
  p = padSize(k);
    
  switch direction
    case 'pre'
       idx{k}   = dimNums(mod(-p:M-1, M) + 1);
    
    case 'post'
      idx{k}   = dimNums(mod(0:M+p-1, M) + 1);
    
    case 'both'
      idx{k}   = dimNums(mod(-p:M+p-1, M) + 1);
  
  end
end
b = a(idx{:});

%%%
%%% SymmetricPad
%%%
function b = SymmetricPad(a, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  dimNums = [1:M M:-1:1];
  p = padSize(k);
    
  switch direction
    case 'pre'
      idx{k}   = dimNums(mod(-p:M-1, 2*M) + 1);
            
    case 'post'
      idx{k}   = dimNums(mod(0:M+p-1, 2*M) + 1);
            
    case 'both'
      idx{k}   = dimNums(mod(-p:M+p-1, 2*M) + 1);
  end
end
b = a(idx{:});

%%%
%%% ReplicatePad
%%%
function b = ReplicatePad(a, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
  M = size(a,k);
  p = padSize(k);
  onesVector = ones(1,p);
    
  switch direction
    case 'pre'
      idx{k}   = [onesVector 1:M];
            
    case 'post'
      idx{k}   = [1:M M*onesVector];
            
    case 'both'
      idx{k}   = [onesVector 1:M M*onesVector];
  end
end
 b = a(idx{:});

%%%
%%% ParseInputs
%%%
function [a, method, padSize, padVal, direction] = ParseInputs(varargin)

% default values
a         = [];
method    = 'constant';
padSize   = [];
padVal    = 0;
direction = 'both';

% checknargin(2,4,nargin,mfilename);

a = varargin{1};

padSize = varargin{2};
% checkinput(padSize, {'double'}, {'real' 'vector' 'nonnan' 'nonnegative' ...
%                    'integer'}, mfilename, 'PADSIZE', 2);

% Preprocess the padding size
if (numel(padSize) < ndims(a))
    padSize           = padSize(:);
    padSize(ndims(a)) = 0;
end

if nargin > 2

    firstStringToProcess = 3;
    
    if ~ischar(varargin{3})
        % Third input must be pad value.
        padVal = varargin{3};
%        checkinput(padVal, {'numeric' 'logical'}, {'scalar'}, ...
%                   mfilename, 'PADVAL', 3);
        
        firstStringToProcess = 4;
        
    end
    
    for k = firstStringToProcess:nargin
        validStrings = {'circular' 'replicate' 'symmetric' 'pre' ...
                        'post' 'both'};
        string = checkstrs(varargin{k}, validStrings, mfilename, ...
                           'METHOD or DIRECTION', k);
        switch string
         case {'circular' 'replicate' 'symmetric'}
          method = string;
          
         case {'pre' 'post' 'both'}
          direction = string;
          
         otherwise
          error('Images:padarray:unexpectedError', '%s', ...
                'Unexpected logic error.')
        end
    end
end
    
% Check the input array type
if strcmp(method,'constant') && ~(isnumeric(a) || islogical(a))
    id = sprintf('Images:%s:badTypeForConstantPadding', mfilename);
    msg1 = sprintf('Function %s expected A (argument 1)',mfilename);
    msg2 = 'to be numeric or logical for constant padding.';
    error(id,'%s\n%s',msg1,msg2);
end

function out = checkstrs(in, valid_strings, function_name, ...
                         variable_name, argument_position)
%CHECKSTRS Check validity of option string.
%   OUT = CHECKSTRS(IN,VALID_STRINGS,FUNCTION_NAME,VARIABLE_NAME, ...
%   ARGUMENT_POSITION) checks the validity of the option string IN.  It
%   returns the matching string in VALID_STRINGS in OUT.  CHECKSTRS looks
%   for a case-insensitive nonambiguous match between IN and the strings
%   in VALID_STRINGS.
%
%   VALID_STRINGS is a cell array containing strings.
%
%   FUNCTION_NAME is a string containing the function name to be used in the
%   formatted error message.
%
%   VARIABLE_NAME is a string containing the documented variable name to be
%   used in the formatted error message.
%
%   ARGUMENT_POSITION is a positive integer indicating which input argument
%   is being checked; it is also used in the formatted error message.

%   Copyright 1993-2003 The MathWorks, Inc.
%   $Revision: 1.3.4.4 $  $Date: 2003/05/03 17:51:45 $

% Except for IN, input arguments are not checked for validity.

% checkinput(in, 'char', 'row', function_name, variable_name, argument_position);

start = regexpi(valid_strings, ['^' in]);
if ~iscell(start)
    start = {start};
end
matches = ~cellfun('isempty',start);
idx = find(matches);

num_matches = length(idx);

if num_matches == 1
    out = valid_strings{idx};

else
    out = substringMatch(valid_strings(idx));
    
    if isempty(out)
        % Convert valid_strings to a single string containing a space-separated list
        % of valid strings.
        list = '';
        for k = 1:length(valid_strings)
            list = [list ', ' valid_strings{k}];
        end
        list(1:2) = [];
        
        msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                       upper(function_name), num2ordinal(argument_position), ...
                       variable_name);
        msg2 = 'to match one of these strings:';
        
        if num_matches == 0
            msg3 = sprintf('The input, ''%s'', did not match any of the valid strings.', in);
            id = sprintf('Images:%s:unrecognizedStringChoice', function_name);
            
        else
            msg3 = sprintf('The input, ''%s'', matched more than one valid string.', in);
            id = sprintf('Images:%s:ambiguousStringChoice', function_name);
        end
        
        error(id,'%s\n%s\n\n  %s\n\n%s', msg1, msg2, list, msg3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = substringMatch(strings)
%   STR = substringMatch(STRINGS) looks at STRINGS (a cell array of
%   strings) to see whether the shortest string is a proper substring of
%   all the other strings.  If it is, then substringMatch returns the
%   shortest string; otherwise, it returns the empty string.

if isempty(strings)
    str = '';
else
    len = cellfun('prodofsize',strings);
    [tmp,sortIdx] = sort(len);
    strings = strings(sortIdx);
    
    start = regexpi(strings(2:end), ['^' strings{1}]);
    if isempty(start) || (iscell(start) && any(cellfun('isempty',start)))
        str = '';
    else
        str = strings{1};
    end
end