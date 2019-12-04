function info = spaFITSInfo(filename)
%
% FITSINFO  Information about FITS file
%
% Produces results compatible to the Matlab fitsinfo routine
% but does not support all the cases, only the file formats
% used by SPARTA. Differently from the Matlab routine, unsigned data
% types are natively supported and some related problems of the
% original Matlab routine are solved here.

vv=ver;
octave=0;
for k=1:size(vv,2)
	if (strcmp(vv(k).Name,'Octave')==1)
	    octave=1;
	end
end	    

if nargin ~= 1
   error('nargin ~= 1');
end

%Open file 
fid = fopen(filename,'r','ieee-be');
if fid==-1
  error('spaFITSInfo:1','Cannot open "%s" for reading.', filename);
else
  info(1).Filename = fopen(fid);
end

if (octave==1)
    d = stat(info.Filename);
    info.FileModDate = d.mtime;
    into.FileSize = d.size;
    info.Contents = {};
else    
    d = dir(info.Filename);
    info.FileModDate = d.date;
    info.FileSize = d.bytes;
    info.Contents = {};
end

%Read primary header.  All FITS files have a primary header 
[info.PrimaryData,headertype,datasize] = readHDU(fid);

%Check for SIMPLE keyword in all standard fits files
if (~strcmp('SIMPLE',info.PrimaryData.Keywords(:,1)))
  error('spaFITSInfo:2','"%s" is not a standard FITS file.',filename);
else
  info.Contents{1} = 'Primary';
end

while feof(fid)==0
% go to the next HDU until EOF is reached
% data must be multiple of 2880 bytes
    if datasize ~= 0
        if rem(datasize,2880) == 0
            bytes_to_seek = datasize;
        else
            bytes_to_seek = 2880 - rem(datasize,2880) + datasize;
        end
        %bytes_to_seek = 2880 - rem(datasize,2880) + datasize;
        status = fseek(fid,bytes_to_seek,0);
        if status == -1
            warning('spaFITSInfo:1','seek failed because of incorrect header information.',filename);
            status = fclose(fid);
            return;
        end
    end
 
    [extensioninfo,headertype,datasize,eof] = readHDU(fid);
    if eof==0
        switch headertype
            case 'binary'
                if (~isfield(info,'BinaryTable'))
                    info.BinaryTable = extensioninfo;
                else
                    info.BinaryTable(end+1) = extensioninfo;
                end
                info.Contents{end+1} = 'Binary Table';
                info.BinaryTable(end).ExtensionOffset = info.BinaryTable(end).Offset+ info.BinaryTable(end).RowSize*info.BinaryTable(end).Rows;
        end % switch
    end % if eof==0
end

%Close file
status = fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% readHDU
%% Read a header unit from a FITS file. File cursor is placed at the end of
%% the header (an integral number of 2880 bytes).

function [info,headertype,datasize,eof] = readHDU(fid)
% avoid warnings
eof=0;
info = struct();
headertype = '';
datasize=0;

% read the first keyword of the HDU and classify it
[card, readcount] = fread(fid,80,'uchar');
if feof(fid)
    eof=1;
    return;
elseif readcount ~= 80
   status = fclose(fid);
   error('FITS card image was not 80 bytes. FITS file may be invalid or corrupt.');
end

card =  [char(card)]';
[keyword,value,comment] = parseFITScard(card);

switch lower(keyword)
 case 'simple'
  headertype = 'primary';
  info = struct('DataType',[],'Size',[],'DataSize',[],'MissingDataValue',[],'Intercept',[],'Slope',[],'Offset',[],'Keywords',{});
 case 'xtension'
  switch lower(value(~isspace(value)))
   case 'bintable'
    headertype = 'binary';
    info = struct('Rows',[],'RowSize',[],'NFields',[],'FieldFormat',{},'FieldPrecision',{},'FieldSize',[],'DataSize',[],'MissingDataValue',{},'Intercept',[],'Slope',[],'Offset',[],'ExtensionSize',[],'ExtensionOffset',[],'Keywords',{});
  end
end

info(1).Keywords{1,1} = keyword;
info(1).Keywords{1,2} = value;
info(1).Keywords{1,3} = comment;

while (feof(fid)==0)
    [card, readcount] = fread(fid,80,'uchar');
    if feof(fid)
        return;
    elseif readcount ~= 80
        status = fclose(fid);
        error('FITS card image was not 80 bytes. FITS file may be invalid or corrupt.');
    end

    card =  [char(card)]';
    [keyword,value,comment] = parseFITScard(card);
    
    switch lower(headertype)
        case 'primary'
            info = processPrimaryKeywords(info,keyword,value);
        case 'binary'
            info = processBinaryKeywords(info,keyword,value);
    end
    
    info.Keywords{end+1,1}=keyword;
    info.Keywords{end,2}=value;
    info.Keywords{end,3}=comment;
    if strcmp(keyword,'END')==1
        break
    end
end

info.DataSize = 0;
switch headertype
 case 'primary'
  if ~isempty(info.Size)
    info.DataSize = prod(info.Size)*derivePrecision(info.DataType)/8;
  end
 case 'binary'	  
  info.DataSize = info.RowSize*info.Rows;
end

%Move file position to integral number of 2880 bytes
curr_pos = ftell(fid);
if rem(curr_pos,2880) == 0
     bytes_to_seek = 0;
else
     bytes_to_seek = 2880 - rem(curr_pos,2880);
end

%datasize needs to be output seperately because of the possible extension data
%for binary data
if strcmp(headertype,'binary')
  datasize = info.DataSize+info.ExtensionSize;
else
  datasize = info.DataSize;
end

status = fseek(fid,bytes_to_seek,0);
info.Offset = ftell(fid);
if status == -1
  error('spaFITSInfo:3','%s','seek failed, file is corrupted');
end

% END readHDU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parseFITScard
% Read keyword, value and comment from a FITS card image
% "The keyword shall be a left justified, 8-character, blank-filled, ASCII
%  string with no embedded blanks."

function [keyword, value, comment] = parseFITScard(card)

keyword = card(1:8);
keyword = keyword(~isspace(keyword));
card = card(9:end);

%Separate Value / Comment
% "If a comment is present, it must be preceded by a slash"...
slashidx = findstr('/',card);

%Find quotes, but first remove double single quotes.
idx = findstr(card,'''''');
tempcard = card;
for i=1:length(idx)
  tempcard(idx(i):idx(i)+1) = '';
end
quoteidx = findstr('''',tempcard);
if isempty(slashidx)
  value = card;
  comment = '';
else
  if length(quoteidx)>1 
    %Value is a quoted char string, might need to fix value / comment
    if (and( (slashidx(1) > quoteidx(1)) , (slashidx(1) < quoteidx(2))))
      %Slash was part of the value, find 1st slash after the last quote
      slashidx = slashidx(slashidx>quoteidx(2));
      if isempty(slashidx)
	%No comment
	value = card;
	comment = '';
      else
	value = card(1:(slashidx(1)-1));
	comment = card(slashidx(1):end);
      end
    else
      value = card(1:(slashidx(1)-1));
      comment = card(slashidx(1):end);
    end
  else
    value = card(1:(slashidx(1)-1));
    comment = card(slashidx(1):end);
  end
end

  if (and( (~isempty(comment)) , (length(comment)>=2)))
  %Get rid of / character, not part of the value or comment.
  %Trailing blanks are not mentioned in the standard so leave them alone
  comment = comment(2:end); 
end

%Value
if strcmp(value(1:2),'= ')
  %Keyword is a value keyword
  value = value(3:end);
  %Remove leading/trailing blanks from value
  % must use sscanf to solve the case of 'value   ' to obtain 'value'
  value = sscanf(value,' %s');
  if strcmp(value(1),'''')
    %Quoted character string. Begin and end quotes are not part of
    %the character string.  Trailing spaces are not significant.
    value = deblank(value);
    if strcmp(value(end),'''')
      value = value(2:(end-1));
    else
      %No closing quote, keep reading though.
    end
    %Replace double single quotes with single quote
    value(findstr(value,'''''')) = '';
  elseif ~isempty(findstr(value(1),'('))
      error('spaFITSInfo:4','%s','complex values not supported');
  elseif ~isempty(str2num(value))
    %Numeric value.  Integers and floating point numbers are treated the same
    value = str2num(value);
elseif (or( (strcmp(value(value~=' '),'T')) , (strcmp(value(value~=' '),'F'))))
    %Logical value. Remove spaces.
    value = value(value~=' ');
  else
    warning('spaFITSInfo:2','%s','Unknown format on header');
  end
else
  %Comment keyword.
  %Return value as ascii text. Remove trailing blanks.
  value = deblank(value);
  return;  
end
%End parseFITScard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processPrimaryKeywords
%
%

function [info,extensions] = processPrimaryKeywords(info,keyword,value)
extensions = 0;
switch lower(keyword)
 case 'bitpix'
  info(1).DataType =  deriveDataType(value);
  info(1) = processUnsigned(info);
 case 'naxis'
 case 'extend'
   if strcmp(value,'T')
     extensions = 1;
   end
 case 'bzero'
  info(1).Intercept = value;
  info(1) = processUnsigned(info);
 case 'bscale'
  info(1).Slope = value;
  info(1) = processUnsigned(info);
 case 'blank'
  info(1).MissingDataValue = value;
 case 'end'
 otherwise
  %Take care of NAXISN keywords
  if findstr('naxis',lower(keyword))
    dim = sscanf(keyword(6:end),'%f');
    info(1).Size(dim) = value;
    if length(info(1).Size)==1
      info(1).Size(2) = 1;
    end
  end    
end    
%END processPrimaryKeywords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deriveDataType
%
%

function dataType = deriveDataType(value)
switch value
 case 8
  dataType = 'int8';
 case 16 
  dataType = 'int16';
 case 32
  dataType = 'int32';
 case -32
  dataType = 'single';
 case -64
  dataType = 'double';
end
%END deriveDataType

function precision = derivePrecision(dataType)
switch dataType
 case 'int8'
  precision = 8;
 case 'uint8'
  precision = 8;
 case 'int16' 
  precision = 16;
 case 'uint16' 
  precision = 16;
 case 'int32'
  precision = 32;
 case 'uint32'
  precision = 32;
 case 'single'
  precision = 32;
 case 'double'
  precision = 64;
end
%END derivePrecision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processUnsigned
%
% FITS does not define unsigned data types, but instead uses BZERO/BSCALE
% Determine whether a data type is unsigned

function info=processUnsigned(info)

switch info(1).DataType
    case 'int8'
        if (and( (info(1).Slope == 1) , ( info(1).Intercept == 128 )))
        info(1).DataType = 'uint8';
        end
    case 'int16'
        if (and( (info(1).Slope == 1) , (info(1).Intercept == 32768)))
        info(1).DataType = 'uint16';
        end
    case 'int32'
        if (and( (info(1).Slope == 1) , (info(1).Intercept == 2147483648)))
        info(1).DataType = 'uint32';
        end
end
% END processUnsigned

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processBinaryKeywords
%
%

function [info,extensions] = processBinaryKeywords(info,keyword,value)
extensions = 0;
endFound = 0;
switch lower(keyword)
 case 'end'
  endFound = 1;
 case 'naxis1'
  info.RowSize = value;
 case 'naxis2'
  info.Rows = value;
 case 'tfields'
  info.NFields = value;
  info.Slope = ones(1,value);
  info.Intercept = zeros(1,value);
  info.MissingDataValue = cell(1,value);
 case 'pcount'
  info.ExtensionSize = value;
 otherwise
   if findstr('tscal',lower(keyword));
    tscale_idx = sscanf(keyword(6:end),'%i');
    info.Slope(tscale_idx) = value;
   elseif findstr('tzero',lower(keyword));
    tzero_idx = sscanf(keyword(6:end),'%i');
    info.Intercept(tzero_idx) = value;
   elseif findstr('tnull',lower(keyword));
     tnull_idx = sscanf(keyword(6:end),'%i');
     info.MissingDataValue{tnull_idx} = value;
   elseif findstr('tform',lower(keyword))
    idx = sscanf(keyword(6:end),'%s');
    repeat = sscanf(value,'%f',1);
    if isempty(repeat)
      repeat = 1;
    end
    format = sscanf(value,' %*i%c',1);
    if isempty(format)
      format = sscanf(value,' %c',1);
    end
    % The value for tformN is a format defined in the FITS standard.  The
    % form is rTa, where r is the repeat, T is the precision and a is a
    % number undefined by the standard. 
    %
    % Binary Table Format valid T (precision) characters and # bytes
    %
    %  L - Logical             1
    %  X - Bit                 *
    %  B - Unsigned Byte       1
    %  I - 16-bit integer      2
    %  J - 32-bit integer      4
    %  A - Character           1
    %  E - Single              4
    %  D - Double              8
    %  C - Complex Single      8
    %  M - Complex Double      16
    %  P - Array Descriptor    8
    switch format
     case 'L'
      format = 'char';
     case 'X'
      format = ['bit' num2str(repeat+8-rem(repeat,8))];
      if repeat~=0
	repeat = 1;
      end
     case 'B'
      format = 'uint8';
     case 'I'
      format = 'int16';
     case 'J'
      format = 'int32';
     case 'A'
      format = 'char';
     case 'E'
      format = 'single';
     case 'D'
      format = 'double';
     case 'C'
      format = 'single complex';
     case 'M'
      format = 'double complex';
     case 'P'
      format = 'int32';
      if repeat~=0
	repeat = 2;
      end
    end
    info.FieldFormat{str2num(idx)} = value;
    info.FieldPrecision{str2num(idx)} = format;
    info.FieldSize(str2num(idx)) = repeat;
  end    
end    
%END KNOWNBINTABLEKEYWORDS
