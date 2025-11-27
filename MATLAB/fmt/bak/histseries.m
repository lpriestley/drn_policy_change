function histseries(tseries, i, j, k, varargin)

hist(tseries,sqrt(length(tseries)));

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 
