function [data,params,hdr] = spaRTDFRead(filename,forceDataType)
% function [data,params,hdr] = spaRTDFRead(filename)
%
% RCS: "@(#) $Id: spaRTDFRead.m,v 0.18 2011/06/03 16:06:31 vltsccm Exp $"
%
% (c) ESO by Enrico Fedrigo, modified by Rob Donaldson

vv=ver;
octave=0;
for k=1:size(vv,2)
	if (strcmp(vv(k).Name,'Octave')==1)
	    octave=1;
            v=vv(k).Version;
            if (v(1:3) == '3.8')
	      octave=0;
            end
	end
end	

hdr=spaFITSInfo(filename);
data={};
params={};

if ~isfield(hdr,'BinaryTable')
    disp(['File "' filename '" is not an RTDF file']);
    return
end
                
numberOfField = hdr.BinaryTable(1).NFields;
data = struct('info','SPARTA Data');
names = {};

for x=1:numberOfField
    keyword=['TTYPE' num2str(x) ];
    
    for y=1:length(hdr.BinaryTable(1).Keywords)
        if (strcmp(hdr.BinaryTable(1).Keywords(y,1),keyword) == 1 )
            [names{x}] = deal(hdr.BinaryTable(1).Keywords(y,2));
        end
    end
end

rows = hdr.BinaryTable(1).Rows;

num = {};
type = {};
offset = {}; % TZERO

localdata = {};
for x=1:length(hdr.BinaryTable(1).FieldFormat)
    ff=char(hdr.BinaryTable(1).FieldFormat(x));
    num{end+1}=str2num(ff(1:end-1));
    type{end+1}=ff(end);
    offset{end+1}=hdr.BinaryTable(1).Intercept(x);
    localdata{x}=zeros(rows,num{x});
end

numdata = length(num);

fid = fopen(filename,'r','b');
fseek(fid,hdr.BinaryTable(1).Offset,-1);
updateStep = fix(rows/100);
if (updateStep>0)
    if octave
        waitbar(0);
    else
        bar = waitbar(0,['Reading file ' filename]);
        pos=get(bar,'Position');
        pos(2)=pos(2)-50;
        pos(3)=275;
        pos(4)=50;
        set(bar,'Position',pos);
    end
end

for n=1:rows
    for x=1:numdata
        switch(type{x})
        case 'E'
            readdata=fread(fid,num{x},'float32');
        case 'J'
            readdata=fread(fid,num{x},'uint32');
        case 'I'
            readdata=fread(fid,num{x},'uint16');
        end
        localdata{x}(n,:)=readdata-offset{x}; % must be modified by TZERO
    end
    
    if (updateStep>0)
        if octave
            waitbar(n/rows);
        else
            if (rem(n,updateStep)==0)
                waitbar(n/rows,bar);    
            end
        end
    end
end

if (updateStep>0)
    if octave
        waitbar(1);
    else
        close(bar);
    end
end
fclose(fid);

%
% Assemble structure
%

for x=1:numberOfField
    if (nargin() > 1)
        if (forceDataType == 1)
            switch(type{x})
            case 'E'
                localdata{x}=single(localdata{x});
            case 'J'
                localdata{x}=uint32(localdata{x});
            case 'I'
                localdata{x}=uint16(localdata{x});
        end
        end
    end
    data=setfield(data,char(names{x}),localdata{x});   
end

% parse primary header
numberOfParameters = length(hdr.PrimaryData.Keywords);

params = {};

for x=1:numberOfParameters
    type=hdr.PrimaryData.Keywords{x,1};

    if (strcmp(type,'HIERARCH') == 1 )
        keyword = hdr.PrimaryData.Keywords{x,2};
        params{end+1}=keyword;
    elseif (strcmp(type,'VERSION') == 1)
        version = hdr.PrimaryData.Keywords{x,2};
    end
end

clear localdata;
