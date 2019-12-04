function [res header]=readbin(fn)
%function [res header]=readbin(fn)
    if istype(fn, '.fits') || istype(fn, '.fits.gz')
        isfits=1;
    elseif istype(fn, '.bin') || istype(fn, '.bin.gz')
        isfits=0;
    end
    if ~exist(fn)
        if exist([fn, '.bin'])
            fn=[fn, '.bin'];
            isfits=0;
        elseif exist([fn, '.fits'])
            fn=[fn, '.fits'];
            isfits=1;
        else
            error(sprintf('%s not found\n', fn));
        end
    end

    %MATLAB does not support transparent gunzip with fopen. 
    %So try to decompress file first
    try
        fn2=gunzip(fn); fn2=fn2{1};
        if ~istype(fn, '.gz')
            movefile(fn2, fn);
            fn2=fn;
        end
        fprintf('Decompressed %s to %s\n', fn, fn2);
        fn=fn2;
    end
    if isfits
        fid=fopen(fn,'rb','b');%handle big endian
    else
        fid=fopen(fn,'rb');
    end
    if isfits %we assembly multiple arrays to a cell
        count=0;
        err=0;
        while ~feof(fid) && ~err
            count=count+1;
            [res{count,1} header{count,1} err]=readbin_do(fid, isfits);
        end
        %last cell is fake when err happens
        if err
            res=res(1:end-1);
            header=header(1:end-1);
        end
        if length(res)==1
            res=res{1};
            header=header{1};
        end
    else
        [res header err]=readbin_do(fid, isfits);
    end
    fclose(fid);
function res=istype(fn, suffix)
    if length(fn)>length(suffix) && strcmp(fn(end-length(suffix)+1:end), suffix)
        res=1;
    else
        res=0;
    end
function [magic, nx, ny, header]=readfits_header(fid)
%Every header line is 80 bytes. Total 36 lines or 2880 bytes.
    END=0;
    page=0;
    bitpix=0;
    nx=0;
    ny=0;
    header='';
    while ~END
        start=1;
        if page==0
            res=fread(fid, 80, 'char*1=>char*1')';
            if length(res)==0 || feof(fid)
                break;
            end
            if ~strcmp(res(1:6), 'SIMPLE') && ~strcmp(res(1:16), 'XTENSION= ''IMAGE')
                disp('Garbage in fits file');
                break;
            end
            res=fread(fid, 80, 'char*1=>char*1')';
            bitpix=str2num(res(11:30));
            res=fread(fid, 80, 'char*1=>char*1')';
            naxis=str2num(res(11:30));
            if(naxis>2)
                error('Data type not supported\n');
            end
            if naxis>0
                res=fread(fid, 80, 'char*1=>char*1')';
                nx=str2num(res(11:30));
            end
            if naxis>1
                res=fread(fid, 80, 'char*1=>char*1')';
                ny=str2num(res(11:30));
            end
            start=start+3+naxis;
        end
        for i=start:36
            res=fread(fid, 80, 'char*1=>char*1')';
            if strcmp(res(1:7), 'COMMENT')
                header=[header res(11:max(find(~isspace(res))))];
            elseif strcmp(res(1:3), 'END')
                END=1;
            end
        end
        page=page+1;
    end
    switch (bitpix)
     case -32
      magic=25608;
     case -64
      magic=25602;
     case 32
      magic=25621;
     case 16
      magic=25611;
     case 8
      magic=25610;
     case 0
      magic=0;
    end

function [magic, nx, ny, header]=readbin_header(fid)
    M_SKIP=26112;
    M_HEADER=25856;

    magic=fread(fid,1,'uint32=>uint32');
    if magic==M_SKIP
        magic=fread(fid,1,'uint32=>uint32');
    end
    header='';
    while magic==M_HEADER
        nlen=fread(fid, 1, 'uint64=>uint64');
        header=[header char(fread(fid, nlen, 'char*1=>char*1')')];
        nlen2=fread(fid, 1, 'uint64=>uint64');
        magic2=fread(fid, 1, 'uint32=>uint32');
        if nlen~=nlen2 || magic2~=M_HEADER
            error('Header verification failed\n');
        end
        magic=fread(fid,1,'uint32=>uint32');
        if magic==M_SKIP
            magic=fread(fid,1,'uint32=>uint32');
        end
    end
    nx=fread(fid,1,'uint64=>double');
    ny=fread(fid,1,'uint64=>double');
function res=fread_fits(varargin)
    res=fread(varargin{:});
    fid=varargin{1};
    tmp=whos('res');
    nbyte=tmp.bytes;
    nbyte2=ceil(nbyte/2880)*2880;
    if nbyte<nbyte2 %read in and discard the padding
        tmp=fread(fid, nbyte2-nbyte, 'char*1=>char*1');
    end
function [res header err]=readbin_do(fid, isfits)
    M_CSP64=25600;
    M_SP64=25601;
    M_DBL=25602;
    M_INT64=25603;
    M_CMP=25604;
    M_INT32=25605;
    M_CSP32=25606;
    M_SP32=25607;
    M_FLT=25608;
    M_ZMP=25609;
    M_INT8=25610;
    M_INT16=25611;
    MC_CSP=25616;
    MC_SP=25617;
    MC_DBL=25618;
    MC_INT64=25619;
    MC_CMP=25620;
    MC_INT32=25621;
    
    MCC_ANY=25633;
    MCC_DBL=25634;
    MCC_CMP=25636;


    MAT_SP=65281;
    MAT_CSP=65282;
    if isfits
        [magic, nx, ny, header]=readfits_header(fid);
        zfread=@fread_fits;
    else
        [magic, nx, ny, header]=readbin_header(fid);
        zfread=@fread;
    end
    if magic<=0
        err=1;
        res=[];
        return;
    else
        err=0;
    end
    if nx==0 || ny==0
        res=[];
        return;
    end
    convert_sparse=0;
    switch magic
     case {MCC_ANY, MCC_DBL, MCC_CMP, MC_CSP, MC_SP, MC_DBL, MC_CMP, MC_INT32, MC_INT64}
      res=cell(nx,ny);
      header2=cell(nx,ny);
      for ii=1:nx*ny
          [res{ii} header2{ii} err]=readbin_do(fid, isfits);
          if err
              break;
          end
      end
      header2{end+1}=header;
      header=header2;
     case {M_SP64, M_CSP64}
      nz=zfread(fid,1,'uint64=>uint64');
      Jc=zfread(fid,ny+1,'uint64=>uint64');
      Ir=zfread(fid,nz,'uint64=>uint64');
      if magic==M_CSP64
          P=zfread(fid,[2 nz],'double=>double');
      else
          P=zfread(fid,nz,'double=>double');
      end
      convert_sparse=1;
     case {M_SP32, M_CSP32}
      nz=zfread(fid,1,'uint64=>uint64');
      Jc=zfread(fid,ny+1,'uint32=>uint32');
      Ir=zfread(fid,nz,'uint32=>uint32');
      if magic==M_CSP32
          P=zfread(fid,[2 nz],'double=>double');
      else
          P=zfread(fid,nz,'double=>double');
      end
      convert_sparse=1;
     case M_DBL
      res=zfread(fid,[nx,ny],'double=>double');
     case M_FLT
      res=zfread(fid,[nx,ny],'float=>float');
     case M_INT64
      res=zfread(fid,[nx,ny],'int64=>int64');
     case M_INT32
      res=zfread(fid,[nx,ny],'int32=>int32');
     case M_CMP
      tmp=zfread(fid,[2*nx,ny],'double=>double');
      res=squeeze(tmp(1:2:end,:)+1i*tmp(2:2:end,:));
      clear tmp; 
     case M_ZMP
      tmp=zfread(fid,[2*nx,ny],'float=>float');
      res=squeeze(tmp(1:2:end,:)+1i*tmp(2:2:end,:));
      clear tmp; 
     otherwise
      fprintf('magic=%d\n',magic);
      error('invalid magic number.');
    end
    if convert_sparse
        Ic=zeros(nz,1);
        for iy=1:ny
            Ic(Jc(iy)+1:Jc(iy+1))=iy;
        end
        Ir=Ir+1;%start from 1 now.
        if size(P,1)==2
            Pr=P(1,:)';
            Pi=P(2,:)';
            P=Pr+1i*Pi;
        end
        res=sparse(Ir,Ic,P,nx,ny);
    end
    
