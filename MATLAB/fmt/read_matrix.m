function [m] = read_matrix(fname)
% [m] = READ_MATRIX(fname)
%
%  Reads a matrix from a file in the binary format written by
%  write_binary_matrix() in miscmaths (FSL library)
%
%  See also: write_matrix (for binary matrices), load (for ascii matrices)

% open file in big-endian
endian='b';
fid=fopen(fname,'r','b');
testval = fread(fid,1,'uint32');
% check if this gives the correct magic number
if (testval~=42),
  fclose(fid);
  % otherwise try little-endian
  fid=fopen(fname,'r','l');
  endian='l';
  testval = fread(fid,1,'uint32');
  if (testval~=42),
    disp('Can not read this file format');
    return;
  end
end

	% ditch the padding
  dummy=fread(fid,1,'uint32');
	% read the number of rows and columns
  nrows=fread(fid,1,'uint32');
  ncols=fread(fid,1,'uint32');
  m=fread(fid,nrows*ncols,'double');
  m=reshape(m,nrows,ncols);
fclose(fid);
return;

