bps = dbstatus;
vars = who;

for i = 1:length(vars)
    if ~( strcmp(vars{i}, 'bps') || strcmp(vars{i}, 'vars') )
        eval(['clear ' vars{i}]);
    end
end
clear i;
clear vars;

dbstop(bps);

close all
clc

pause(0.1);