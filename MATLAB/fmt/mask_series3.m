function out = mask_series3(data,mask)
% function out = mask_series3(data,mask)
%
% Takes a 4D dataset and returns the 2D matrix of the timeseries for the
% mask (same as vols2matrix, but faster for a small ROI)

[xVals yVals zVals] = find3(mask);
out = zeros(length(xVals),size(data,4));

for i=1:length(xVals)
    out(i,:) = squeeze(data(xVals(i),yVals(i),zVals(i),:));
end
