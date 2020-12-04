%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSGCSF_FANUC class file for manipulators 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef PSGCFS_FANUC
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
       MAX_I_ITER = 1;
       epsilon_I = 1e-4;
       
       % calculation statistics
       iter_I = 1;
       iter_O = 1;
       total_iter = 0;
       
   end
   
   methods
       %%% _init_
       function self = PSGCFS_FANUC(val,val2,varargin)
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
                % reset
                self.eval.u_old = self.u;               
                self.iter_I = 1;
                % get feasible set
                self = self.get_con();
                % do PSG
                self = self.inner_PSG_5(); 
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
       function self = inner_PSG_5(self)
          % do PSG
          while self.stop_inner()~=true
              self.eval.cost_old = self.eval.cost_new;
              self = self.PSG_update_arm();
              % get cost
              self.eval.cost_new = self.eval.get_cost(self.u);
              self.iter_I = self.iter_I +1;
          end
          % get new xref
          self.x_=[];
          for i=2:self.sys_info.H+1
              self.sys_info.xR(:,i)=self.sys_info.robot.A([1:self.sys_info.njoint,self.sys_info.njoint+2:1+2*self.sys_info.njoint],[1:self.sys_info.njoint,self.sys_info.njoint+2:1+2*self.sys_info.njoint])*self.sys_info.xR(:,i-1)...
                  +self.sys_info.robot.B([1:self.sys_info.njoint,self.sys_info.njoint+2:1+2*self.sys_info.njoint],1:self.sys_info.nu)*self.u((i-2)*self.sys_info.nu+1:(i-1)*self.sys_info.nu);
              self.x_=[self.x_;self.sys_info.xR(:,i)];
          end
          self.eval.x_ = self.x_;
       end
       
       %%% One PSG step  
       function self = PSG_update_arm(self)
           u_ = self.u;
           % Stochastic gradient descent
           u_ = u_ - self.sys_info.alpha*(self.dcostArm_f() + 10*normrnd(0,0.1,[self.nn,1])./((self.iter_O)^2 + 1));
           % project to feasible set
           self = self.Projection(u_);         
       end
       
       %%% Project SGD to feasible set 
       function self = Projection(self,u_)
           exitflag = 2;
           H = eye(self.nn);
           f = -u_;
           %[u_proj,cost,exitflag,output,lambda] = quadprog(H,f,self.Ainq,self.binq,[],[],-0.2*ones(self.nn,1),0.2*ones(self.nn,1),[],self.options);
           [u_proj,cost,exitflag,output,lambda] = quadprog(H,f,self.Ainq,self.binq,[],[],[],[],[],self.options);
           
           self.total_iter = self.total_iter + output.iterations;
                   %    while exitflag~=1
                %        self = 
                %        [u_proj,fval,exitflag,output] = quadprog(H,f,self.Ainq,self.binq,[],[],-0.2*ones(self.nn,1),0.2*ones(self.nn,1),u_,options);
                %    end
           self.u  = u_proj;
       end
       
       %%% Gradient direction from the cost function
       function df = dcostArm_f(self)
           df = self.sys_info.QQ*self.u + self.sys_info.ff;
       end
      
       %%% Inner-loop stopping condition
       function [STOP_I] = stop_inner(self)
           STOP_I  = false;           
           delta = norm(self.eval.cost_new - self.eval.cost_old);
           if delta < self.epsilon_I || self.iter_I > self.MAX_I_ITER
               STOP_I  = true;
           end                    
       end
        
       %%% get feasible set
       function self = get_con(self)
           DH = self.sys_info.robot.DH;
           base = self.sys_info.robot.base;
           nstate = self.sys_info.nstate;
           njoint = self.sys_info.njoint;          

           Lstack=[];Sstack=[];
           for j = 1:size(self.obs,2)
               I=[];       
               f = @(x) self.dist_arm_all(x,base,self.obs{j}.l,self.sys_info.robot);
               for i=1:self.sys_info.H
                   theta=self.x_(nstate*(i-1)+1:nstate*(i-1)+njoint);
                   [distance,linkid] = self.dist_arm_all(theta,base,self.obs{j}.l,self.sys_info.robot);
                   I = [I;distance-self.obs{j}.D];
%                    Diff=zeros(njoint,1);
%                    for s=1:njoint
%                        [Diff(s),~]=derivest(@(x) self.dist_arm_id([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,self.obs{j}.l,self.sys_info.robot.cap,linkid),theta(s),'Vectorized','no');
%                    end
                   Diff = num_jac(f,theta); Diff = Diff';                   
                   Bj=self.sys_info.Baug((i-1)*nstate+1:i*nstate,1:self.nn);
                   s=I(i)-Diff'*Bj(1:njoint,:)*self.u;
                   l=-Diff'*Bj(1:njoint,:);
                   
%                    Bj=Baug((i-1)*nstate+1:i*nstate,1:horizon*nu);
%                    s=I(i)-Diff'*Bj(1:njoint,:)*uref;
%                    l=-Diff'*Bj(1:njoint,:);
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