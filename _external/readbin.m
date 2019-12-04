function [res header]=readbin(fn0)
    fn=fn0;
    if length(fn)>6 &&strcmp(fn(end-6:end),'.bin.gz')%end with bin.gz
        if ~exist(fn,'file')
            fn=fn(1:end-3);
        end
    elseif length(fn)>3 &&strcmp(fn(end-3:end),'.bin') %end with bin
        if ~exist(fn,'file')
            fn=[fn '.gz'];
        end
    else %no suffix
        fn=[fn '.bin'];
        if ~exist(fn,'file')
            fn=[fn '.gz'];
        end
    end
    if ~exist(fn)
        error(sprintf('%s not found\n', fn));
    end
    
    if strcmp(fn(end-2:end),'.gz')
        fprintf('uncompressing %s\n',fn);
        system(sprintf('gunzip -f %s',fn));
        fn=fn(1:end-3);
    end
    fid=fopen(fn,'rb');
    [res header]=readbin_do(fid);
function [res header]=readbin_do(fid, header)
    M_CSP64=25600;
    M_SP64=25601;
    M_CSP32=25606;
    M_SP32=25607;
    M_DBL=25602;
    M_INT64=25603;
    M_CMP=25604;
    M_INT32=25605;
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
    M_HEADER=25856;
    M_SKIP=26112;
    magic=fread(fid,1,'uint32');
    if magic==M_SKIP
        magic=fread(fid,1,'uint32');
    end
    header='';
    while magic==M_HEADER
        nlen=fread(fid, 1, 'uint64');
        header=[header char(fread(fid, nlen, 'char*1')')];
        nlen2=fread(fid, 1, 'uint64');
        magic2=fread(fid, 1, 'uint32');
        if nlen~=nlen2 || magic2~=M_HEADER
            error('Header verification failed\n');
        end
        magic=fread(fid,1,'uint32');
        if magic==M_SKIP
            magic=fread(fid,1,'uint32');
        end
    end
    nx=fread(fid,1,'uint64');
    ny=fread(fid,1,'uint64');
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
          [res{ii} header2{ii}]=readbin_do(fid);
      end
      header2{end+1}=header;
      header=header2;
     case {M_SP64, M_CSP64}
      nz=fread(fid,1,'uint64');
      Jc=fread(fid,ny+1,'uint64');
      Ir=fread(fid,nz,'uint64');
      if magic==M_CSP64
          P=fread(fid,[2 nz],'double');
      else
          P=fread(fid,nz,'double');
      end
      convert_sparse=1;
     case {M_SP32, M_CSP32}
      nz=fread(fid,1,'uint64');
      Jc=fread(fid,ny+1,'uint32');
      Ir=fread(fid,nz,'uint32');
      if magic==M_CSP32
          P=fread(fid,[2 nz],'double');
      else
          P=fread(fid,nz,'double');
      end
      convert_sparse=1;
     case M_DBL
      res=fread(fid,[nx,ny],'double');
     case M_INT64
      res=fread(fid,[nx,ny],'int64');
     case M_INT32
      res=fread(fid,[nx,ny],'int32');
     case M_CMP
      tmp=fread(fid,[2*nx,ny],'double');
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
    