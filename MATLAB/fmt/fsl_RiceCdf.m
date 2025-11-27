function p = fsl_RiceCdf(x,A,s2)
% p = fsl_RiceCdf(ordinates,A,sigma2)
%
% Returns the cumulative probabilities for a 
% (set of) ordinates given a Rician distribution 
% with amplitude A and noise parameter s2.
%
% Uses Simpson quadrature to evaluate the integral
% of the Rician pdf. I suspect it could be done
% in a more elegant way.
%

if any(x)<0, error('Rician distribution only defined for positive values'); end
if A<0, error('Rician distribution only defined for positive amplitudes'); end

p = zeros(size(x));
for i=1:length(x)
  p(i) = quadl(@(y) fsl_RicePdf(y,A,s2),0,x(i));
end

return