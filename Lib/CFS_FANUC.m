%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CSF_FANUC class file for manipulators 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef CFS_FANUC
   properties
       % setup parameters
       obs cell
       sys_info struct
       nn {mustBeNumeric}
       ROBOT = 'M16iB'                          % 'M16iB' or 'M200i'
      
       % outputs
       u {mustBeNumeric}                   % velocity       
       x_ {mustBeNumeric}                  % trajectory
                   
       % constraint
       Ainq {mustBeNumeric}
       binq {mustBeNumeric}
       
       % functions
       dist_arm_all = @(th,base,obs,robot) dist_arm_3D_Heu_2(th,base,obs,robot);
       
       % sover settings       
       options = optimoptions('quadprog','Display','off');
       
       % evaluation setting
       eval EVAL
       
       % calculation statistics       
       iter_O = 1;
       total_iter = 0; 
       
   end
   
   methods
       %%% _init_
       function self = CFS_FANUC(val,val2,varargin)
            % get problem info
            self.obs = val;
            self.sys_info = val2;                
            self.nn = self.sys_info.H*self.sys_info.nu;
            if ~isempty(varargin)
                self.ROBOT = varargin{1};                    
            end 
            % initialization
            switch( self.ROBOT)
                case 'M200i'
                    self.dist_arm_all = @(th,base,obs,robot) dist_arm_3D_200i_2(th,base,obs,robot);
                case '2L'
                    self.dist_arm_all = @(th,base,obs,robot) dist_arm_2L(th,base,obs,robot);                    
            end       
            self.x_ = self.sys_info.x_;
            self.u = zeros(self.nn,1); 
            % initialize evaluation functions
            self.eval = EVAL(val2);
       end
       
       %%% main function
       function self = optimizer(self)          
           self.eval.cost_new = self.eval.get_cost(self.u);
           while self.eval.stop_outer(self.iter_O)~=true
               % save
               self.eval.u_old = self.u;    
               self.eval.cost_old = self.eval.cost_new;             
               % get convex feasible set
               self = self.get_con();
               % solve QP
               self = self.Solve_QP();
               % get cost
               self.eval.cost_new = self.eval.get_cost(self.u);
               % store results
               self.eval = self.eval.store_result(self.u);
               % next step                
               self.iter_O = self.iter_O+1;
           end
       end  
       
       %%%%%%%%%%%%%%%% supporting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% PSG inner-loop with  
       function self = Solve_QP(self)
          % solve QP
          [unew ,cost,exitflag,output,lambda]= quadprog(self.sys_info.QQ,self.sys_info.ff,self.Ainq,self.binq,[],[],-self.sys_info.MAX_input,self.sys_info.MAX_input,[],self.options);
          self.u = unew;          
          % get new xref
          self.eval.x_old = self.x_;
          self.x_=[];
          for i=2:self.sys_info.H+1
              self.sys_info.xR(:,i)=self.sys_info.robot.A([1:self.sys_info.njoint,self.sys_info.njoint+2:1+2*self.sys_info.njoint],[1:self.sys_info.njoint,self.sys_info.njoint+2:1+2*self.sys_info.njoint])*self.sys_info.xR(:,i-1)...
                  +self.sys_info.robot.B([1:self.sys_info.njoint,self.sys_info.njoint+2:1+2*self.sys_info.njoint],1:self.sys_info.nu)*self.u((i-2)*self.sys_info.nu+1:(i-1)*self.sys_info.nu);
              self.x_=[self.x_;self.sys_info.xR(:,i)];
          end
          self.eval.x_ = self.x_;
          % get total iteration
          self.total_iter = self.total_iter + output.iterations;
       end
       
       %%% get feasible set
       function self = get_con(self)
           DH = self.sys_info.robot.DH;
           base = self.sys_info.robot.base;
           nstate = self.sys_info.nstate;
           njoint = self.sys_info.njoint;           
           %%% probably won't use
           theta_all = reshape(self.x_,[nstate,self.sys_info.H]);
           [cap_all] = get_cap_all(theta_all(1:njoint,:),base,self.sys_info.robot);
           Lstack=[];Sstack=[];
           for j = 1:size(self.obs,2)
               I=[];
               f = @(x) self.dist_arm_all(x,base,self.obs{j}.l,self.sys_info.robot);
               for i=1:self.sys_info.H
                   theta=self.x_(nstate*(i-1)+1:nstate*(i-1)+njoint);
                   [distance,linkid] = self.dist_arm_all(theta,base,self.obs{j}.l,self.sys_info.robot);
                   
                   I = [I;distance-self.obs{j}.epsilon];
                   Diff = num_jac(f,theta); Diff = Diff';                   
                   Bj=self.sys_info.Baug((i-1)*nstate+1:i*nstate,:);
                   s=I(i)-Diff'*Bj(1:njoint,:)*self.u;
                   l=-Diff'*Bj(1:njoint,:);                  

                   Sstack=[Sstack;s];
                   Lstack=[Lstack;l];
                   % velosity constraints
                   Lstack=[Lstack;self.sys_info.Baug((i-1)*nstate+njoint+1:i*nstate,:)];
                   Sstack=[Sstack;self.sys_info.lim-self.sys_info.Aaug((i-1)*nstate+njoint+1:i*nstate,:)*self.sys_info.xR(:,1)];
                   Lstack=[Lstack;-self.sys_info.Baug((i-1)*nstate+njoint+1:i*nstate,:)];
                   Sstack=[Sstack;self.sys_info.lim+self.sys_info.Aaug((i-1)*nstate+njoint+1:i*nstate,:)*self.sys_info.xR(:,1)];
                 
               end
           end
           self.Ainq = Lstack;
           self.binq = Sstack;
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   end
end