function [obs, sim] = run_swat(x, Extra)

path = Extra.path_to_python;
id =  char(java.util.UUID.randomUUID);

x= sprintf('%.3f,', x);
x = x(1:end-1); % strip final comma

command = sprintf('%s %s ''%s'' ''%s''', path, Extra.wrapper, x, id);
system(command);

outfile = [id,'.dat'];
output = load(outfile);
obs = output(:,1); sim = output(:,2);

delete(outfile)