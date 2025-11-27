function cock
% See  woolrich
list={'Behrens','Woolrich','Woolrich','Woolrich','Beckmann','Smith','Flitney','Johansen-Berg','Mortimer'};
disp('Thinking');
for i=1:3,
      disp('.');
      pause(0.5);
end;
display(char(list(floor(rand*length(list))+1)));


if(0)
[s,w]=dos('whoami');
if(strcmp(w(1:end-1),'woolrich'))
  disp('Thinking');
  pause(1)
display('Behrens');
elseif(strcmp(w(1:end-1),'behrens'))
   disp('Thinking');
  pause(1)
display('Woolrich');
else
   disp('Thinking');
  pause(1)
display('Steve Smith');
end;
end;




