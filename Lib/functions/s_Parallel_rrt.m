%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel RRT (script function)
% With planning algorithm: RRT, RRT*
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


self_ = {};
self = [];
num_seed = 6;                      % depends on #cores of the computer 
path_fail = true(num_seed,1);
iter_rrt = 0;
while all(path_fail)
    routeL = 1000*ones(num_seed,1);
    parfor i=1:num_seed 
    self_{i} = RRT_FANUC(obs,sys_info,goalxyz,region_g,region_s,sample_off,'M200i','RRT');
    self_{i} = self_{i}.find_route();
    path_fail(i) = self_{i}.fail;
        if self_{i}.fail ~=true
            routeL(i) = size(self_{i}.route,2);    
        end
    end
    iter_rrt = iter_rrt+1;
end
%%
[path_length,id] = min(routeL);
self = self_{id};


