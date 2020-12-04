%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOMP_FANUC class file for manipulators 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef CHOMP_FANUC
   properties
       % setup parameters
       obs cell
       sys_info struct
       nn {mustBeNumeric}
       ROBOT = 'M16iB'                          % 'M16iB' or 'M200i'
      
       % outputs
       u {mustBeNumeric}                   % velocity
       x_ {mustBeNumeric}                  % trajectory
       
       % functions
       dist_arm_all = @(th,DH,base,obs,cap) dist_arm_3D_Heu(th,DH,base,obs,cap);
       dist_arm_id = @(th,DH,base,obs,cap,id) dist_link_Heu(th,DH,base,obs,cap,id);
       
       % evaluation setting
       eval EVAL
       
       % calculation statistics       
       iter_O = 1;
       total_iter = 0;       
       
   end
   
   methods
       %%% _init_
       function self = CHOMP_FANUC(val,val2,uu,varargin)
            % get problem info
            self.obs = val;
            self.sys_info = val2;                
            self.nn = self.sys_info.H*self.sys_info.nu;
            if ~isempty(varargin)
                self.ROBOT = varargin{1};                    
            end 
            % initialization
            if self.ROBOT == 'M200i'
                self.dist_arm_all = @(th,DH,base,obs,cap) dist_arm_3D_200i(th,DH,base,obs,cap);
                self.dist_arm_id = @(th,DH,base,obs,cap,id) dist_link_200i(th,DH,base,obs,cap,id);
            end       
            self.x_ = self.sys_info.x_;            
            self.u = uu;
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
               % do CHOMP               
               self = self.CHOMP_update_arm();
               % get cost
               self.eval.cost_new = self.eval.get_cost(self.u) + self.fobs_m();
               % store results
               self.eval = self.eval.store_result(self.u);
               % next step                 
               self.iter_O = self.iter_O+1; 
           end
       end  
       
       %%%%%%%%%%%%%%%% supporting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% PSG inner-loop with  
       function self = CHOMP_update_arm(self)
          u_ = self.u;
          self.u = u_ - self.sys_info.alpha*3*(self.dcostArm_f() + 2000*self.dcostObs_f());   
          % new ref
          self.x_=[];
          for i=2:self.sys_info.H+1
              self.sys_info.xR(:,i)=self.sys_info.robot.A([1:self.sys_info.njoint,7:6+self.sys_info.njoint],[1:self.sys_info.njoint,7:6+self.sys_info.njoint])*self.sys_info.xR(:,i-1)...
                  +self.sys_info.robot.B([1:self.sys_info.njoint,7:6+self.sys_info.njoint],1:self.sys_info.nu)*self.u((i-2)*self.sys_info.nu+1:(i-1)*self.sys_info.nu);
              self.x_=[self.x_;self.sys_info.xR(:,i)];
          end
       end
       
       %%% Gradient direction from the cost function
       function df = dcostArm_f(self)
           df = self.sys_info.QQ*self.u + self.sys_info.ff;
       end
       
       %%% get cost 
       function c_all = fobs_m(self)
           nstate = self.sys_info.nstate;
           njoint = self.sys_info.njoint;  
           c_all = 0;
           for i=1:self.sys_info.H
               theta=self.x_(nstate*(i-1)+1:nstate*(i-1)+njoint);    
               for j = 1:self.obs{1}.num_obs
               Dfx = self.dm_f(theta,j);               
                   for s =1:njoint
                       if Dfx(s)<0
                           c_x = -Dfx(s)+(1/2)*self.obs{j+1}.epsilon;    
                       elseif Dfx(s)>=0 && Dfx(s)<=self.obs{j+1}.epsilon
                           c_x = (1/(2*self.obs{j+1}.epsilon))*(Dfx(s)-self.obs{j+1}.epsilon)^2;
                       else
                           c_x = 0;
                       end
                       c_all = c_all + c_x;
                   end
               end
           end

       end
       
       %%% distance function
       function [d] = dm_f(self,theta,j)
           DH = self.sys_info.robot.DH;
           base = self.sys_info.robot.base;
           njoint = size(theta,1);
           for i=1:njoint
               DH(i,1)=theta(i);
           end%theta,d,a,alpha
           d = [];
           if size(base,2)>1
               base=base';
           end
           pos=CapPos(base,DH,self.sys_info.robot.cap);           
           for i=1:njoint
               [dis, points] = distLinSeg(pos{i}.p(:,1),pos{i}.p(:,2), self.obs{j+1}.l(:,1),self.obs{j+1}.l(:,2));
               if norm(dis)<0.0001
                   dis = -norm(points(:,1)-pos{i}.p(:,2));
               end        
               d(i) =dis - self.obs{j+1}.D;
           end           
       end
       
       %%% get obstacle cost gradient
       function dc_all = dcostObs_f(self)
           DH = self.sys_info.robot.DH;
           base = self.sys_info.robot.base;
           nstate = self.sys_info.nstate;
           njoint = self.sys_info.njoint;
           dc_all = zeros(self.nn,1);
           for i=1:self.sys_info.H
               theta=self.x_(nstate*(i-1)+1:nstate*(i-1)+njoint);    
               for j = 1:self.obs{1}.num_obs
                   Dfx = self.dm_f(theta,j);
                   [dis, linkid] = min(Dfx);
                   dDfx = zeros(njoint,1);
                   if Dfx(linkid)<0
                       for s=1:njoint
                           [dDfx(s),~]=derivest(@(x) self.dist_arm_id([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,self.obs{j+1}.l,self.sys_info.robot.cap,linkid),theta(s),'Vectorized','no');
                       end            
                       dc_x = -[dDfx'*self.sys_info.Baug((i-1)*njoint+1:i*njoint,:)]';
                   elseif Dfx(linkid)>=0 && Dfx(linkid)<=self.obs{j+1}.epsilon
                       for s=1:njoint
                           [dDfx(s),~]=derivest(@(x) self.dist_arm_id([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,self.obs{j+1}.l,self.sys_info.robot.cap,linkid),theta(s),'Vectorized','no');
                       end 
                       dc_x = [(1/self.obs{j+1}.epsilon)*(Dfx(linkid)-self.obs{j+1}.epsilon)*dDfx'*self.sys_info.Baug((i-1)*njoint+1:i*njoint,:)]';
                   else
                       dc_x =  zeros(self.nn,1);
                   end
                   dc_all  = dc_all + dc_x;           
               end
           end
        end
       
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   end
end