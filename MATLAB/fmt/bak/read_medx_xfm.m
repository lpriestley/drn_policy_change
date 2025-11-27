function [g,a,s] = read_medx_xfm(fname)
% [g,a,s] = READ_MEDX_XFM(fname)
%
%  Extracts the GenericReslice matrix (g),
%   the AlignLinearReslice matrix (a), and
%   the ShadowTransform matrix (s)  (or UserTransform) from a
%   MEDx transformation file (given by fname)
%  Any 4x4 or 3x4 matrices should be formatted
%  Other, unrecognised, sizes will be an nx1 vector
%
g=[];
a=[];
s=[];
fid=fopen(fname);
if (fid<0),
  error(['File ',fname,' not found.']);
end
while (~feof(fid)),
  str=fscanf(fid,'%s',1);
  if (strcmp(str,'/GenericReslice')),
    tstg=get_matrix(fid);
    if length(tstg>0),  g=tstg;  end
  end
  if (strcmp(str,'/MotionCorrectionReslice')),
    tsta=get_matrix(fid);
    if length(tsta>0),  a=tsta;  end
  end
  if (strcmp(str,'/ShadowTransform')),
    tsts=get_matrix(fid);
    if length(tsts>0),  s=tsts;  end
  end
end
fclose(fid);
return;


%-----------------------------------------

function g = get_matrix(fid)

g=[];
str1 = fscanf(fid,'%s',2);
if (~strcmp(str1,'<</matrix')),
  return;
end
str1 = fscanf(fid,'%s',1);
if (~strcmp(str1,'[')),
 disp('Error in reading Reslice Matrix');
end
 str1 = fscanf(fid,'%s',1);
 while (~strcmp(str1,']')),
   g=[g str2num(str1)];
%   disp(['READ ',str1]);
   str1 = fscanf(fid,'%s',1);
 end
if (length(g)==16),
  g=reshape(g,4,4).';
elseif (length(g)==12),
  g=reshape(g,4,3).';
end
return;
