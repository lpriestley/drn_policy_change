function linewhisker(mn,err,colour)
if(nargin==2)
  colour='r';
end

ho;
ind=1:length(mn);
line([ind-0.2; ind+0.2],[mn; mn],'LineWidth',3,'Color',colour);
line([ind; ind],[mn-err; mn+err],'LineWidth',1,'Color','k')
line([ind-0.1; ind+0.1],[mn+err; mn+err],'LineWidth',1,'Color', 'k')
line([ind-0.1; ind+0.1],[mn-err; mn-err],'LineWidth',1,'Color','k');