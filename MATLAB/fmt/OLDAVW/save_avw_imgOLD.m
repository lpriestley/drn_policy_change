function save_avw_img(img,fname,vtype);
%  SAVE_AVW_IMG(img,fname,vtype)
%
%  Save an array (img) as an analyse file (only the .img) 
%   for either a 2D or 3D or 4D array (automatically determined)
%  Note: 
%        assumes the data uses a MEDx coordinate convention
% 
%  vtype is a single character string: 'b' byte, 's' short or 'f' float
%
%  See also: SAVE_AVW, SAVE_AVW_HDR, READ_AVW, READ_AVW_HDR, READ_AVW_IMG
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
fnimg=strcat(fname,'.img');

fp=fopen(fnimg,'w');
dims = size(img);

%% DEFUNCT
%% flip y dimension to be consistent with MEDx
%% dat=flipdim(img,2);

dat = img;
dat = reshape(dat,prod(dims),1);

switch vtype
  case 'f'
    vtype2='float';
  case 's'
    vtype2='short';
  case 'b'
    vtype2='schar';
  case 'u'
    vtype2='uchar';
end;

fwrite(fp,dat,vtype2);
fclose(fp);

