%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PGDCFS on manipulator 
% Robot model: FANUC M200i
%
% Jessica Leu+k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

%% Load robot info
ROBOT = 'M200i';
robot = robotproperty2(ROBOT);
njoint=5;nstate=10;nu=5;DH=robot.DH;
offset= [3150, 8500,330]./1000;  % center of the scene
robot.base = offset';

%% Initialization  
load('data/200i_xori.mat') 
x0 = route_wp(:,1);                 % get starting th's 
xg = route_wp(:,end);               % get end th's
xR=[];
xR(:,1)=[x0; zeros(5,1)];
horizon = 30;                    % specify horizon/ time duration
% Generatie reference
% line in the state space
x_ = [];
for i=1:nstate/2
    x_ = [x_; linspace(x0(i),xg(i),horizon+1)];
end%
W = zeros(5,horizon+1);
x_ = [x_;W];

xref_ = [];
for i = 1:horizon+1%
    xref_ = [xref_; x_(:,i)];
end   
xori=xref_(nstate+1:end);
x_ = xori;
% zeros input
uref = zeros(horizon*njoint,1);
uori = zeros(horizon*njoint,1);

%% Environment
% Set up obstacles
obs_{1}.num_obs = 1;
obs=[[3606;8313;1]./1000 [3606;8313;1038]./1000];%
D=0.2;
obs_{2}.l = obs;
obs_{2}.D = D;

%% Cost Fn Parameters
Aaug=[];Baug=zeros(horizon*nstate,horizon*nu);Qaug=zeros(horizon*nstate);
Q=[];
Q(1:njoint,1:njoint)=[10 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];
Q(1:njoint,njoint+1:2*njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,1:njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,njoint+1:2*njoint)=[10 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];
for i=1:horizon
    Aaug=[Aaug;robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])^i];
    Qaug((i-1)*nstate+1:i*nstate,(i-1)*nstate+1:i*nstate)=Q*0.1;
    if i==horizon
        Qaug((i-1)*nstate+1:i*nstate,(i-1)*nstate+1:i*nstate)=Q*10000;
    end
    for j=1:i
        Baug((i-1)*nstate+1:i*nstate,(j-1)*nu+1:j*nu)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])^(i-j)*robot.B([1:njoint,7:6+njoint],1:nu);
    end
end
R=eye(horizon*nu);
for i=1:horizon
    R((i-1)*nu+1:i*nu,(i-1)*nu+1:i*nu)=[5 0 0 0 0;
        0 4 1 0 0;
        0 1 2 0 0;
        0 0 0 2 0;
        0 0 0 0 1];
end
R=R+R';
QQ = Baug'*Qaug*Baug+R.*0.1;
ff = ((Aaug*xR(:,1)-xori)'*Qaug*Baug)';

self.QQ = QQ;
self.ff = ff;
self.sys_info.QQ = QQ;
self.sys_info.ff = ff;
self.sys_info.robot = robot;
self.sys_info.H = horizon;
self.sys_info.nstate = nstate;
self.sys_info.njoint = njoint;
self.sys_info.xR = xR;
self.sys_info.nu = nu;


sys_info.Aaug = Aaug;
sys_info.Baug = Baug;
sys_info.QQ = QQ;
sys_info.ff = ff;
sys_info.robot = robot;
sys_info.H = horizon;
sys_info.nstate = nstate;
sys_info.njoint = njoint;
sys_info.xR = xR;
sys_info.nu = nu;
sys_info.x_ = x_;
sys_info.alpha = 1/max(svd(QQ)); 
% solver
sys_info.epsilon_O = 1e-6;

self = PSGCFS_FANUC(obs_,sys_info,'M200i');
self = self.optimizer();
uref = self.u;
x_ = self.x_;


% self = [];
% 
% self.QQ = QQ;
% self.ff = ff;
% self.sys_info.QQ = QQ;
% self.sys_info.ff = ff;
% self.sys_info.robot = robot;
% self.sys_info.H = horizon;
% self.sys_info.nstate = nstate;
% self.sys_info.njoint = njoint;
% self.sys_info.xR = xR;
% self.sys_info.nu = nu;
% self.sys_info.Aaug = Aaug;
% self.sys_info.Baug = Baug;
% self.x_ = x_;
% self.obs = obs_;
% self.dist_arm_all = @(th,DH,base,obs,cap) dist_arm_3D_200i(th,DH,base,obs,cap);
% self.dist_arm_id = @(th,DH,base,obs,cap,id) dist_link_200i(th,DH,base,obs,cap,id);
%                 
% 
% 
% % Parameter
% self.nn = horizon*nu;
% % Solver setup
% self.sys_info.alpha = 1/max(svd(QQ));    % stochastic ratio 
% self.MAX_O_ITER = 10;
% self.MAX_I_ITER = 20;
% self.epsilon = 1e-6; 
% self.options = optimoptions('quadprog','Display','off');
% % Initialize
% self.iter_O = 0;
% self.total_iter = 0;
% self.u = uref;
% self.u_old = ones(self.nn,1);
% self.cost_old = 1000000;
% self.cost_new = get_cost_arm(self);
% cost_all = [];
% 
% tic
% while stop_outer(self)~=true
%     % reset
%     self.u_old = self.u;    
%     self.cost_new = get_cost_arm(self);
%     self.iter_I = 1;
%    
%     % get feasible sert
%     %[self.Ainq,self.binq,Link] = get_con(self.sys_info.robot.DH,robot.base,obs,robot,self.x_,self.u,horizon,nstate,njoint,D,nu,Baug);
%     %[self.Ainq,self.binq,Link] = get_con(self.sys_info.robot.DH,self.sys_info.robot.base,self.obs{2}.l,self.sys_info.robot,self.x_,self.u,self.sys_info.H,self.sys_info.nstate,self.sys_info.njoint,self.obs{2}.D,self.sys_info.nu,self.sys_info.Baug);
%     self = get_con(self);
%     % do PSG
%     self = inner_PSG_5(self);
%     
%     % get cost
%     cost1 = self.u'*self.sys_info.QQ*self.u + self.sys_info.ff'*self.u;
%     cost_all = [cost_all cost1];
%     
%     % get new xref
%     self.x_=[];
%                   for i=2:self.sys_info.H+1
%                       self.sys_info.xR(:,i)=self.sys_info.robot.A([1:self.sys_info.njoint,7:6+self.sys_info.njoint],[1:self.sys_info.njoint,7:6+self.sys_info.njoint])*self.sys_info.xR(:,i-1)...
%                           +self.sys_info.robot.B([1:self.sys_info.njoint,7:6+self.sys_info.njoint],1:self.sys_info.nu)*self.u((i-2)*self.sys_info.nu+1:(i-1)*self.sys_info.nu);
%                       self.x_=[self.x_;self.sys_info.xR(:,i)];
%                   end    
%        
%     
%     % next step
%     self.iter_O = self.iter_O+1;
% end 
% toc
% uref = self.u;
% x_ = self.x_;
%% Visualization 

% plot inputs
uref_ = reshape(uref,[5,horizon]);
uori_ = reshape(uori,[5,horizon]);
figure(2);
plot(uref_','--')
legend('\alpha_1','\alpha_2','\alpha_3','\alpha_4','\alpha_5')
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Input [rad/s^2]')

% plot positions
x_ = reshape(x_,[10,horizon]);
figure(3)
plot(x_(1:5,:)')
hold on
legend('1','2','3','4','5')
% plot ref
xori_ = reshape(xori,[10,horizon]);
plot(xori_(1:5,:)','--')
legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5', ...
        '\theta_{1,ref}','\theta_{2,ref}','\theta_{3,ref}','\theta_{4,ref}','\theta_{5,ref}')
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Input [rad]')

% robot animation
xref_ = reshape(x_,[10,horizon]);
figure(1)
DrawMap2;
delete(h)
[all_ee] = plot_FANUC_200i(obs,xref_,robot,njoint,D,horizon,offset);

%% plot cost all
figure(4)
plot(cost_all)
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Cost')

%% functions

function self = inner_PSG_5(self)           
   while stop_inner(self)~=true
       self.cost_old = self.cost_new;
       self = PSG_update_arm(self);
   end
end


function self = PSG_update_arm(self)
   u_ = self.u;
   u_ = u_ - self.sys_info.alpha*(dcostArm_f(self) + 1*normrnd(0,0.1,[self.nn,1])./((self.iter_I)^2 + 1));
               
self = Projection(self,u_);  
   self.cost_new = get_cost_arm(self);
   self.iter_I = self.iter_I +1;

end


function self = Projection(self,u_)
          exitflag = 2;
          H = eye(self.nn);
          f = -u_;   
          [u_proj,cost,exitflag,output,lambda] = quadprog(H,f,self.Ainq,self.binq,[],[],-0.2*ones(self.nn,1),0.2*ones(self.nn,1),[],self.options);
           self.total_iter = self.total_iter + output.iterations;
           %    while exitflag~=1
        %        self = 
        %        [u_proj,fval,exitflag,output] = quadprog(H,f,self.Ainq,self.binq,[],[],-0.2*ones(self.nn,1),0.2*ones(self.nn,1),u_,options);
        %    end
          self.u  = u_proj;
end  


   function df = dcostArm_f(self)
       df = self.sys_info.QQ*self.u + self.sys_info.ff;
   end

   function cost = get_cost_arm(self)
      cost = self.u'*self.sys_info.QQ*self.u + self.sys_info.ff'*self.u;

   end 
% 

function [STOP_I] = stop_inner(self)
   STOP_I  = false;           
   delta = norm(self.cost_new - self.cost_old);
   if delta < self.epsilon || self.iter_I > self.MAX_I_ITER
       STOP_I  = true;
   end                    
end
       
function [STOP_O] = stop_outer(self)
   STOP_O  = false;           
   if self.iter_O > self.MAX_O_ITER
       STOP_O  = true;
   end                    
end

function self = get_con(self)
           DH = self.sys_info.robot.DH;
           base = self.sys_info.robot.base;
           nstate = self.sys_info.nstate;
           njoint = self.sys_info.njoint;          
           
           Lstack=[];Sstack=[];
           for j = 1:self.obs{1}.num_obs
               I=[];               
               for i=1:self.sys_info.H
                   theta=self.x_(nstate*(i-1)+1:nstate*(i-1)+njoint);
                   [distance,linkid] = self.dist_arm_all(theta,DH(1:njoint,:),base,self.obs{j+1}.l,self.sys_info.robot.cap);
                   I = [I;distance-self.obs{j+1}.D];
                   Diff=zeros(njoint,1);
                   for s=1:njoint
                       [Diff(s),~]=derivest(@(x) self.dist_arm_id([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,self.obs{j+1}.l,self.sys_info.robot.cap,linkid),theta(s),'Vectorized','no');
                   end
                   Bj=self.sys_info.Baug((i-1)*nstate+1:i*nstate,1:self.nn);
                   s=I(i)-Diff'*Bj(1:njoint,:)*self.u;
                   l=-Diff'*Bj(1:njoint,:);        
                   Sstack=[Sstack;s];
                   Lstack=[Lstack;l];
               end
           end
           self.Ainq = Lstack;
           self.binq = Sstack;
       end

