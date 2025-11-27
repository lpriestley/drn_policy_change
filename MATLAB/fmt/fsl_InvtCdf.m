function x = fsl_InvtCdf(p,df)
% x = fsl_InvtCdf(cumulative probabilities,df)
%
% Returns the inverse cumulative probabilities for a 
% (set of) ordinates given a t-distribution 
% with df degrees of freedom.
%
% Solves the inverse problem numerically using
% Newton-Raphson. I suspect there are more direct
% approaches, but this was easy ;-)
%

maxiter = 100;
tol = 1e-6;

x = zeros(size(p));
for i=1:length(p)
   xo = 0;
   for j=2:maxiter
      xn = xo - (fsl_tCdf(xo,df) - p(i)) / fsl_tPdf(xo,df);
      if abs(xn-xo) < tol 
         break;
      else
         xo = xn;
      end
   end
   x(i) = xn;
end

return;