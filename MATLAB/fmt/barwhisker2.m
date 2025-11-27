function barwhisker2(mn,err,colour,ind,width)
if(nargin<3) colour='k';end
if(nargin<4) ind=1:length(mn); end
if(nargin<5) width=0.8; end
bar(ind, mn,width,colour);
ho;

line([ind; ind],[mn; mn+err],'LineWidth',3,'Color','k');
line([ind-0.1; ind+0.1],[mn+err; mn+err],'LineWidth',3,'Color','k');