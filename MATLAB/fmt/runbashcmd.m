function runbashcmd(cmd)

global homedir;

tmp2=sprintf('%s\\matlab\\fmt\\runbash "%s"', homedir,cmd);
dos(tmp2);