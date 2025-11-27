function p = fsl_NPdf(x,mu,s2)
% p = fsl_NPdf(ordinates,mean,variance)
%
% Returns the probabilities for a (set of)
% ordinates given a normal distribution with
% mean mu and variance s2.
%

if s2<=0 , error('Normal-distribution only defined for positive variance'); end

p = (1/sqrt(2*pi*s2)) * exp(-((x-mu).^2)/(2*s2));

return