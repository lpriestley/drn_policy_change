function x = thresh(x, t)

% x = thresh(x, t)
% set all x>t to t and all x<t to 0

x2 = squash(x);
x2(find(x<t)) = 0;
x = reshape(x2,size(x));