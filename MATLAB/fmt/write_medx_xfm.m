function write_medx_xfm(fname,no_voxels,g,a,s,outputdim)
% WRITE_MEDX_XFM(fname,no_voxels,g,a,s,outputdim)
%
%  Writes a MEDx transformation file using
%   the GenericReslice matrix (g),
%   the AlignLinearReslice matrix (a), 
%   the ShadowTransform matrix (s).
%   into a file specified by fname.
%  Any 4x4 or 3x4 matrices will be correctly formatted
%  The outputdim matrix (4x4) gives the sizes for the voxel dimensions
%   e.g. outputdim = diag([0.93 0.93 5 1])
%  The no_voxels gives the volume size in voxels [no_x,no_y,no_z]
%   e.g. no_voxels = [256 256 30]
%
fid=fopen(fname,'w');
fprintf(fid,'%s\n','%!VEST-Transformations');
fprintf(fid,'%s\n','<<');
if (length(g)>0),
     fprintf(fid,'%s\n','    /GenericReslice     <<');
     put_matrix(fid,'matrix',g);
     fprintf(fid,'%s\n','        /order  1');
     put_matrix(fid,'outputsize',no_voxels);
     if (length(outputdim)>0),
      put_matrix(fid,'outputusermatrix',outputdim);
     end
     fprintf(fid,'%s\n','    >>');
end
if (length(a)>0),
     fprintf(fid,'%s\n','    /MotionCorrectionReslice     <<');
     put_matrix(fid,'matrix',a);
     fprintf(fid,'%s\n','        /order  1');
     put_matrix(fid,'outputsize',no_voxels);
     if (length(outputdim)>0),
      put_matrix(fid,'outputusermatrix',outputdim);
     end
     fprintf(fid,'%s\n','        /outputunits    (mm)');
     fprintf(fid,'%s\n','        /talairachcalibrated?   false');
     fprintf(fid,'%s\n','    >>');
end
if (length(s)>0),
     fprintf(fid,'%s\n','    /ShadowTransform     <<');
     put_matrix(fid,'matrix',s);
     fprintf(fid,'%s\n','        /order  1');
     put_matrix(fid,'outputsize',no_voxels);
     if (length(outputdim)>0),
      put_matrix(fid,'outputusermatrix',outputdim);
     end
     fprintf(fid,'%s\n','    >>');
end
fprintf(fid,'%s\n','>>');
fclose(fid);
return;


%-----------------------------------------

function put_matrix(fid,name,g)

if (length(g)<=0),
  return;
end
if (length(g)~=prod(size(g))),
  g=reshape(g.',prod(size(g)),1);
end
fprintf(fid,'%s\n',['        /',name,' [']);
for idx=1:length(g),
  fprintf(fid,'%s\n',['            ',num2str(g(idx))]);
end
fprintf(fid,'%s\n','        ]');
return;
