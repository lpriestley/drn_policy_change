function p = fsl_Chi2Pdf(x,df)
% p = fsl_Chi2Pdf(ordinates,df)
%
% Returns the probabilities for a (set of)
% ordinates given a chi-2 distribution with
% df degrees of freedom.
%

if any(x)<0, error('Chi-2 distribution only defined for positive values'); end
if df<1 || abs(round(df)-df)>eps , error('Chi-2 distribution only defined for df=1,2,...,Inf'); end

p = (1/(2^(df/2)*gamma(df/2))) .* x.^(df/2-1) .* exp(-x/2);

return