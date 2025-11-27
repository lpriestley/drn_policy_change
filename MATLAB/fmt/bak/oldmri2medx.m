function [transmat] = mri2medx(tr_mat,sr,init_size,final_size,finalcentre)
% [transmat] = mri2medx(tr_mat,sr,init_size,final_size,finalcentre)
%
% Converts a transformation matrix between mri format
%   and medx format  (pre and post flipped in y and adjusted for centre)
% It must be given: 
%       resampling matrix (sr)  - a diagonal matrix
%       the sizes of the initial and final
%          volumes init_size, final_size both as [xmax,ymax,zmax].
%             Convention is that xmax = no x voxels , ymax = etc.
%       centre of the final volume (finalcentre) - a 3x1 matrix

csze=size(finalcentre);
if (csze(1)<csze(2)),
  finalcentre=finalcentre(1,1:3).';
end
% posttr transforms from world -> minc voxels of the av305
%xoffset=[final_size(1)-1-finalcentre(1); finalcentre(2:3)];
%posttr=[diag([-1 1 1]), xoffset; 0 0 0 1];
posttr=[diag([1 1 1]), finalcentre; 0 0 0 1];
% pretr transforms from minc voxels -> world
pretr=tr_mat*sr;

flipy=diag([1 -1 1 1]);
% flip1 transforms from medx voxels -> minc voxels (initial volume)
flip1=flipy;
flip1(2,4)=init_size(2)-1;
% flip2 transforms from minc voxels -> medx voxels (final volume)
%flip2=diag([-1 -1 1 1]);
%flip2(1,4)=final_size(1)-1;
flip2=diag([1 -1 1 1]);
flip2(2,4)=final_size(2)-1;

transmat=flip2*posttr*pretr*flip1;
