function [x,y,z] = find3(matrix)
% finds the 3D-cartesian coordinates of the non-zero elements of matrix

if length(size(matrix)) > 3
    error('the given matrix must be 3-D')
end

if size(matrix,3)==1
    [x,y] = find(matrix);
    z = ones(size(x));
else

    xnew = []; ynew = []; znew = [];
    for m = 1:length(matrix(1,1,:))
        [xtemp, ytemp] = find(matrix(:,:,m));
        xnew = [xnew; xtemp];
        ynew = [ynew; ytemp];
        znew = [znew; ones(length(xtemp),1)*m];
    end

    x=xnew; y=ynew; z=znew;
end