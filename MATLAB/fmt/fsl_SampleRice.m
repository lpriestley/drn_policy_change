function s = fsl_SampleRice(n,A,s2)
% sample = fsl_SampleRice(n,A,sigma2)
%
% Returns an nx1 column vector with a
% sample drawn from a Rician distribution 
% with amplitude A and noise parameter s2.
% Uses rand, so resetting rand resets
% fsl_SampleRice.

if A<0, error('Rician distribution only defined for positive or zero amplitudes'); end
if s2<0, error('Rician distribution only defined for positive or zero variances'); end

s = sqrt((A+sqrt(s2)*randn(n,1)).^2 + (sqrt(s2)*randn(n,1)).^2);

return

%
% Silly way of doing things
%
%  p = rand(n,1);
%  s = fsl_InvRiceCdf(p,A,s2);

