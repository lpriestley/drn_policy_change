function plot_mm(mix,data,mm)

% plot_mm(mix,data)
% 
% plots Gaussian Mixture Model compared with data

nmixes = length(mix.centres);
res = (max(data) - min(data))/(length(data)/20);

mind = res*floor(min(data)/res);
maxd = res*ceil(max(data)/res);
al = [mind:res:maxd]';
alneg = [mind:res:-res]';
alplus = [res:res:maxd]';

mixtures = zeros(nmixes,length(al));

hs = hist(data,al);

bar(al,hs/sum(hs));

hold on;

for m=1:nmixes,
  if mix.priors(m) > 0,
    if(mm == 2 & m==3)      
      me=mix.centres(m);
      va=mix.covars(m);
      plus = mix.priors(m)*gammapdf(me^2/va, me/va, alplus)';
      mixtures(m,length(al)-length(alplus)+1:length(al)) = plus; 
      
    elseif(mm == 2 & m==1)
      me=-mix.centres(m);
      va=mix.covars(m);      
      neg = mix.priors(m)*gammapdf(me^2/va, me/va, -alneg)';
      mixtures(m,1:length(alneg)) = neg;      
      
    elseif(mix.covars(m)>0)
      mixtures(m,:) = mix.priors(m)*normpdf(al,mix.centres(m),sqrt(mix.covars(m)))';    
    end;
  end;
end;


for m=1:nmixes,
   
  %plot(al',mix.priors(m)*mixtures(m,:)./sum(sum(mixtures,2)),'r','LineWidth',2);
  plot(al',mixtures(m,:)./sum(sum(mixtures,2)),'r','LineWidth',2);
end;

plot(al',sum(mixtures,1)./sum(sum(mixtures,2)),'y--','LineWidth',2);

hold off;
