%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RRT_FANUC class file for manipulators 
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef RRT_FANUC
   properties
       % setup parameters
       obs cell
       sys_info struct
       pos cell
       goal {mustBeNumeric}
       region_g {mustBeNumeric}
       region_s {mustBeNumeric}
       sample_off {mustBeNumeric}
       dt = 0.2;
       ROBOT = 'M16iB'                          % 'M16iB' or 'M200i'
       SOLVER = 'RRT*'                          % 'RTT*' or 'RRT'
       
       % intermediate states
       parent {mustBeNumeric}
       newNode {mustBeNumeric}       
       isGoalReached logical       
       toNode_dis {mustBeNumeric}               % vector
       
       % outputs
       all_nodes {mustBeNumeric}                % all nodes
       total_dis {mustBeNumeric}
       route {mustBeNumeric}
       fail = 0;
       
       % store data 
       all_ee {mustBeNumeric}                   % all nodes            
           
       % sover settings
       MAX_ITER = 400;  %1000
       bi = 0.5;
       prevLalpha = 0;
              
       % calculation statistics
       node_num = 1;
   end
   
   
   methods
       %%% _init_
       function self = RRT_FANUC(val,val2,val3,val4,val5,val6,varargin)             
                % get problem info
                self.obs = val;
                self.sys_info = val2;                
                self.goal = val3;
                self.region_g = val4;
                self.region_s = val5;
                self.sample_off = val6;
                if ~isempty(varargin)
                    self.ROBOT = varargin{1};
                    self.SOLVER = varargin{2};
                end                   
       end       
       
       %%% Main function       
       function self = find_route(self)
           % initialize tree
           self.newNode = [self.sys_info.x0];
           self.all_nodes = [-1; self.newNode];     
           self.total_dis = [0];
           self = self.goal_reached();
           % generate rrt
           switch self.SOLVER
               case 'RRT'
                   while self.isGoalReached ~=true
                       self = self.getNode();
                       self = self.addNode();               
                       self = self.goal_reached();               
                   end
               case 'RRT*'
                   while self.isGoalReached ~=true
                       self = self.getNode();
                       self = self.addNode();
                       self = self.arrangeNode();               
                       self = self.goal_reached();               
                   end
           end
           % find route
           self.route = self.newNode;
           while self.parent ~= -1
               self.route = [ self.all_nodes(2:end,self.parent) self.route];
               self.parent = self.all_nodes(1,self.parent);
           end
       end      
       
       %%%%%%%%%%%%%%%% supporting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% get feasible node
       function  self = getNode(self)
           self = self.getRandNode();
           [isFeasible, self]= self.feasible();
           while isFeasible ~=true
               self = self.getRandNode();
               [isFeasible, self]= self.feasible();               
           end    
           
       end
       
       %%%% get random node
       function self = getRandNode(self)
           self.toNode_dis = [];
           pp = rand;
           if pp<self.bi
               % specified sample range 
               sampleNode = (rand(self.sys_info.nstate,1)-0.5).*self.region_s*2+self.sample_off;
           else
               sampleNode = self.sys_info.goal_th;
           end

           self.parent = 1;
           dis = norm((self.all_nodes(2:end,1)-sampleNode).*self.sys_info.ratial);
           self.toNode_dis = [dis];
           for i = 2:self.node_num
               
               dis_ = norm((self.all_nodes(2:end,i)-sampleNode).*self.sys_info.ratial);               
               self.toNode_dis = [self.toNode_dis dis_];
               if dis_<dis
                   dis = dis_;
                   self.parent = i;
               end
           end
           
           self.newNode = self.all_nodes(2:end,self.parent)+(sampleNode-self.all_nodes(2:end,self.parent))*0.1/norm(self.all_nodes(2:end,self.parent)-sampleNode);
             
       end
       
       %%% rearrange parent-children pares in RRT*
       function self = arrangeNode(self)
           ids = find(self.toNode_dis <0.2);
           for i = 1:length(ids)
               if self.total_dis(ids(i))>(self.total_dis(end)+self.toNode_dis(ids(i)))
                   self.all_nodes(1,ids(i)) = self.node_num; 
                   self.total_dis(ids(i)) = (self.total_dis(end)+self.toNode_dis(ids(i)));
               end
           end
       end
       
       %%%%%%%%%%%%%%%%%%% supplementary functions %%%%%%%%%%%%%%%%%%%%%%%%
       %%% check random node feasibility
       function [isFeasible, self]= feasible(self)           
            DH = self.sys_info.DH;
            base = self.sys_info.base;            
            isFeasible = true;
            
            for j = 1:size(self.obs,2)
                % load check target
                theta=self.newNode;
                % rearrange 
                for i=1:self.sys_info.nstate
                DH(i,1)=theta(i);
                end 
                if self.ROBOT == 'M200i'
                    DH(2,1) = DH(2,1)- pi/2;
                end
                if size(base,2)>1
                    base=base';
                end
                % get capsule positions
                pos_=CapPos(base,DH,self.sys_info.robot.cap);
                % check distance
                for i=1:self.sys_info.nstate
                    [dis, points] = distLinSeg(pos_{i}.p(:,1),pos_{i}.p(:,2), self.obs{j}.l(:,1),self.obs{j}.l(:,2));
                    if norm(dis)<0.0001
                        dis = -norm(points(:,1)-pos_{i}.p(:,2));
                    end        
                    if dis < self.obs{j}.D
                        isFeasible = false;
                        break;
                    end
                end
                % save coonfiguration
                self.pos = pos_;
                
            end
       end
       
       %%% add feasible mnode to tree
       function self = addNode(self)           
           self.all_nodes = [self.all_nodes [self.parent; self.newNode]];
           self.all_ee = [self.all_ee  self.pos{self.sys_info.nstate}.p(:,1)];
           self.total_dis = [self.total_dis (self.total_dis(self.parent) + self.toNode_dis(self.parent))];
           self.node_num = self.node_num + 1;
           
       end
       
       %%% check if goal region is reached or if MAXiter is reached
       function self = goal_reached(self)           
           self.isGoalReached = false;           
           if (self.goal-self.region_g)<self.newNode
               if self.newNode<(self.goal+self.region_g)
                   self.isGoalReached = true;
               end
           end
           
           if self.node_num>self. MAX_ITER
               disp('Failed to find path.')
               self.fail = true;
               self.isGoalReached = true;
           end
                      
       end
   end
end