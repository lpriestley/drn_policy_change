function save_avw(img,fname,vtype,vsize)
% SAVE_AVW(img,fname,vtype,vsize) 
%
%  Create and save an analyse header (.hdr) and image (.img) file
%   for either a 2D or 3D or 4D array (automatically determined).
%  Note: 
%        assumes the data uses a MEDx coordinate convention
%	 and that slice direction is transverse
%   
%  vtype is a single character string: 'b' byte, 's' short or 'f' float
%  vsize is a vector [x y z tr] containing the voxel sizes in mm and
%  the tr in seconds  (defaults: [1 1 1 3])
%
%  See also: SAVE_AVW_HDR, SAVE_AVW_IMG, READ_AVW, READ_AVW_HDR, READ_AVW_IMG
%

   save_avw_hdr(img,fname,vtype,vsize);
   save_avw_img(img,fname,vtype);
