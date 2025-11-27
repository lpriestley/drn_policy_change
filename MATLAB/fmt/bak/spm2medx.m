function [transmat] = spm2medx(tr_mat,init_size,final_size)
% [transmat] = spm2medx(tr_mat,init_size,final_size)
% Converts a transformation matrix between spm format
%   and medx format  (pre and post flip in x and y)
% It must be given the sizes of the initial and final
%   volumes init_size, final_size both as [xmax,ymax,zmax].
%   Convention is that xmax = no x voxels , ymax = etc.

flip1=[-1 0 0 init_size(1); 0 -1 0 init_size(2); 0 0 1 1; 0 0 0 1];
flip2=[-1 0 0 final_size(1); 0 -1 0 final_size(2); 0 0 1 -1; 0 0 0 1];
transmat=flip2*tr_mat*flip1;
