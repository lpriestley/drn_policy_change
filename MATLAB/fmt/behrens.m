function behrens

   list={'cock','member','huge cock','nob','dick','penis'};

   disp('Thinking');
   for i=1:3,
      disp('.');
      pause(0.5);
   end;
   
   display(char(list(floor(rand*length(list))+1)));