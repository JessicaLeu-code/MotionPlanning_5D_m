%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVAL class file for xy-planning 
%(evaluation functions) 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef EVAL
   properties
       % setup parameters
       obs      cell
       sys_info struct
       nn       {mustBeNumeric}       
      
       % buffer
       u_old {mustBeNumeric}                   % velocity of previous solution
       x_
       x_old {mustBeNumeric}                   % trajectory of previous solution 
            
       % baseline cost
       Cost_b = 0;      
       
       % sover settings       
       MAX_O_ITER = 20;  %2D: 100;
       epsilon_O = 1e-2; %2D: 1e-4;       
       
       % calculation statistics       
       total_iter = 0;
       cost_old = 100000;
       cost_new = 0;
       
       % store results
       cost_all = [];                      % store all cost
       e_cost_all = [];                    % store all cost-change
       e_u_all = [];
       
   end
   methods
       % _init_
       function self = EVAL(val)           
                % get problem info
                self.sys_info = val;                
                self.epsilon_O = self.sys_info.epsilon_O;
                self.MAX_O_ITER = self.sys_info.MAX_O_ITER;
                % initialization
                self.x_ = self.sys_info.x_;
                self.x_old = ones(size(self.x_));
       end
       
       %%%%%%%%%%%%%%%% evaluation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function cost = get_cost(self,u)           
           cost = 0.5*u'*self.sys_info.Qaug*u + self.sys_info.paug'*u+self.sys_info.caug;
       end
       
       function  self = store_result(self,u)
           self.cost_all = [self.cost_all self.cost_new];
           self.e_cost_all = [self.e_cost_all norm(self.cost_old - self.cost_new)];
           self.e_u_all = [self.e_u_all norm(self.u_old - u)];
       end
                   
       function [STOP_O] = stop_outer(self, iter_O)
           STOP_O  = false; 
          %delta = norm(self.cost_new - self.cost_old);
            delta = norm(self.x_ - self.x_old); 
           if delta < self.epsilon_O
               disp(strcat('Converged at step',num2str(iter_O)));
               STOP_O  = true;
           end
           if iter_O>self.MAX_O_ITER
               disp('MAX_ITER');
               STOP_O  = true;
           end
       end
       
       function  cost = get_Cost_b(self)
           [u_b]= quadprog(self.sys_info.Qaug,self.sys_info.paug,[],[],[],[],[],[],[],self.sys_info.options);
           cost = self.get_cost(u_b);           
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
           
       
   end
end