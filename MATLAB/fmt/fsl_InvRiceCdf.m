function x = fsl_InvRiceCdf(p,A,s2)
% x = fsl_InvRiceCdf(cumulative probabilities,A,sigma2)
%
% Returns the inverse cumulative probabilities for a 
% (set of) ordinates given a Rician distribution 
% with amplitude A and noise parameter s2.
%
% Solves the inverse problem numerically using
% Newton-Raphson. 
%

maxiter = 100;
tol = 1e-6;

x = zeros(size(p));
for i=1:length(p)
   xo = A;
   for j=2:maxiter
      xn = xo - (fsl_RiceCdf(xo,A,s2) - p(i)) / fsl_RicePdf(xo,A,s2);
      if abs(xn-xo) < tol 
         break;
      else
         xo = xn;
      end
   end
   x(i) = xn;
end

return;