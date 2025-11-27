function p = fsl_tPdf(x,df)
% p = fsl_tPdf(ordinates,df)
%
% Returns the probabilities for a (set of)
% ordinates given a t distribution with
% df degrees of freedom.
%

if df<1 || abs(round(df)-df)>eps , error('t-distribution only defined for df=1,2,...,Inf'); end

p = (gamma((df+1)/2) / (sqrt(pi*df)*gamma(df/2))) * (1 + x.^2/df).^-((df+1)/2);

return