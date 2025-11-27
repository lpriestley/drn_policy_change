function runcmd(cmd)

global homedir;

if(findstr(homedir,':')>0),
    % windows
    runbashcmd(cmd);
else
    dos(cmd);
end;