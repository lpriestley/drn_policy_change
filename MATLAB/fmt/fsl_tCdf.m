function p = fsl_tCdf(t,df)
% p = fsl_tCdf(ordinates,degrees_of_freedom)
%
% Returns the cumulative probabilities for a 
% (set of) ordinates given a t-distribution 
% with df degrees of freedom.
%

if df<1 || abs(round(df)-df)>eps , error('t-distribution only defined for df=1,2,...,Inf'); end

x = (t+sqrt(t.^2+df)) ./ (2*sqrt(t.^2+df));
p = betainc(x,df/2,df/2);

return