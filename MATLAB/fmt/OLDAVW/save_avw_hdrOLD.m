function save_avw_hdr(img,fname,vtype,vsize)
% SAVE_AVW_HDR(img,fname,vtype,vsize) 
%
%  Create and save an analyse header file
%   for either a 2D or 3D or 4D array (automatically determined).
%  Note: 
%        assumes the data uses a MEDx coordinate convention
%	 and that slice direction is transverse
%   
%  vtype is a single character string: 'b' byte, 's' short or 'f' float
%  vsize is a vector [x y z tr] containing the voxel sizes in mm and
%  the tr in seconds  (defaults: [1 1 1 3])
%
%  See also: SAVE_AVW, SAVE_AVW_IMG, READ_AVW, READ_AVW_HDR, READ_AVW_IMG 
%

% swap first and second argument in case save_avw_img convention is
% used
check=length(size(fname));
if(check~=2)
   tmp=img;
   img=fname;
   fname=tmp;
end

% remove extension if it exists
if ( (length(findstr(fname,'.hdr'))>0) | ...
        (length(findstr(fname,'.img')>0)) ),
  fname=fname(1:(length(fname)-4));
end

% remove input file:
!touch hdrcreate.txt;
!rm hdrcreate.txt;

% remove headerfile
fname2=strcat(fname,'.hdr');
tmpstr1=sprintf('!touch %s', fname2);
tmpstr2=sprintf('!rm %s', fname2);

eval(tmpstr1);
eval(tmpstr2);

% establish dynamic range
imgmax=ceil(max(max(max(max(img)))));
imgmin=floor(min(min(min(min(img)))));

% create file to use as input into header program
dims = [size(img) 1 1];

if(nargin==2)
  vtype='s';
  vsize=[1 1 1 3];
elseif(nargin==3)
  tmp=size(vtype);
  if(tmp(2)==1)
     vsize=[1 1 1 3];
  else
     vsize=vtype;
     if size(vsize,2)==3
	vsize=[vsize 3];
     end;
     vtype='s';
  end
else
  tmp=size(vtype);
  if(tmp(2)==3)
     tmp2=vtype;
     vtype=vsize;
     vsize=tmp2;
  end
end

if (length(vsize)<3),
  vsize(3)=1;
end
if (length(vsize)<4),
  vsize(4)=3;
end

 fid = fopen('hdrcreate.txt','w');
  if imgmin~=imgmax
    fprintf(fid,'%d\n%d\n%d\n%d\n%s\nv\n%6.4f\n%6.4f\n%6.4f\n%6.4f\nr\n%6.0f\n%6.0f\ns\n',dims(1),dims(2), dims(3),dims(4),vtype,vsize(1),vsize(2),vsize(3),vsize(4),imgmin,imgmax);
  else
    fprintf(fid,'%d\n%d\n%d\n%d\n%s\nv\n%6.4f\n%6.4f\n%6.4f\n %6.4f\ns\n',dims(1),dims(2), dims(3),dims(4),vtype,vsize(1),vsize(2),vsize(3), vsize(4));  
  end;
  fclose(fid);

% call header program

tmp=sprintf('! /usr/local/fsl/bin/header -n %s < hdrcreate.txt \n',fname);
eval(tmp);
disp(' ');

% remove input file:
!rm hdrcreate.txt;

