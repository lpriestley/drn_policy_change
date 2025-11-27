function z = t_to_z(t, dof);

% z = t_to_z(t, dof);
   
%for i=1:length(t),
%  z(i) = sign(t(i))*(2^0.5)*erfinv(1-2*tdist(abs(t(i)),dof));
%end;

z = sign(t).*(2^0.5).*erfinv(1-2.*tdist(abs(t),dof));