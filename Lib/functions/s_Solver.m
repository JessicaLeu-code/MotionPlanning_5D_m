%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver (script function)
% With planning algorithm: CHOMP,CFS, PSGCFS
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(SOLVER,2)
    % select solver
    if ~isempty(SOLVER{i})
        switch(SOLVER{i}.name)
            case 'CHOMP'
                % CHOMP
                tic
                self = CHOMP_FANUC(obs,sys_info,uref,ROBOT);
                self = self.optimizer();
                Time_{1,cc} = toc; 
                Cost{1,cc} = self.eval.cost_all;
                e_u_all{1,cc} = self.eval.e_u_all;
                e_cost_all{1,cc} = self.eval.e_cost_all;
            case 'PSGCFS'
                % PSGCFS
                tic
                self = PSGCFS_FANUC(obs,sys_info,ROBOT);
                self = self.optimizer();
                Time_{2,cc} = toc; 
                Cost{2,cc} = self.eval.cost_all;
                e_u_all{2,cc} = self.eval.e_u_all;
                e_cost_all{2,cc} = self.eval.e_cost_all;
            case 'CFS'
                % CFS
                tic
                self = CFS_FANUC(obs,sys_info,ROBOT);
                self = self.optimizer();
                Time_{3,cc} = toc; 
                Cost{3,cc} = self.eval.cost_all;
                e_u_all{3,cc} = self.eval.e_u_all;
                e_cost_all{3,cc} = self.eval.e_cost_all;
                cost_all = self.eval.cost_all;
        end

        % get data
        uref = self.u;
        x_ = self.x_;
        xref_ = reshape(x_,[nstate,horizon]);

        % robot animation
        if SOLVER{i}.plot == true        
            figure(1)
            [all_ee] =plot_2L(obs,xref_,robot,njoint,view_area);
        end
    end
end



