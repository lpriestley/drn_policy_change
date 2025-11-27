function plotseries_reg(tseries, i, j, k, x)
tseries=tseries-mean(tseries);
%tseries=detrend(tseries,'linear',[200;400;600]);
tseries=detrend(tseries,'linear',[1:40:950]);
plot(tseries);




if(nargin>4)
  ho;  
  x=x-mean(x);
%  x=detrend(x,'linear',[200;400;600]);
    x=detrend(x,'linear',[1:40:950]);
  bet=pinv(x)*tseries;
  plot(x*bet,'r');

end;

title(sprintf('Slice %d; X: %d  Y: %d',k,i,j)); 

hold off;