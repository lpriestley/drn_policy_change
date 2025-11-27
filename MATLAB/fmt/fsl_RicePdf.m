function p = fsl_RicePdf(x,A,s2)
% p = fsl_RicePdf(ordinates,A,sigma2)
%
% Returns the probabilities for a (set of)
% ordinates given a Rician distribution with
% amplitude A and noise parameter s2.
%

if any(x)<0, error('Rician distribution only defined for positive values'); end
if A<0, error('Rician distribution only defined for positive amplitudes'); end

% The expression below is equivalent to
% p = x/s2 * exp(-(x.^2+A^2)/(2*s2)) * besseli(0,x*A/s2);
% but a lot less prone to overflow.

logp = log(x) - log(s2) - (A^2 + x.^2)/(2*s2) + log(besseli(0,x*A/s2,1)) + x*A/s2;
p = exp(logp);


return