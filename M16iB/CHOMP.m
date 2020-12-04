%% main CHOMP
clear all 
close all
%% build cost

robot=robotproperty(3);
ROBOT = 'M16iB';                   % 'M200i' or 'M16iB'
robot = robotproperty2(ROBOT);
njoint=5;nstate=10;nu=5;DH=robot.DH;
base=[3250, 8500,0]./1000;

load('data/M16_ref_2.mat')
load('data/good_xori.mat')
x_ = xuori;
xR=[];xR(:,1)=x_(1:nstate);
horizon=size(uref,1)/njoint;
xori=x_(nstate+1:end);
x_=xori;
uori = uref;

obs=[[4106;8313;1]./1000 [4106;8313;1338]./1000];%

D=0.2;
epsilon =  0.05;

obs_{1}.num_obs = 1;
obs_{2}.D = 0.2;
obs_{2}.l = obs;
obs_{2}.epsilon =  0.05;
%% plot ref
xori_ = reshape(xori,[10,horizon]);
figure;plot(xori(1:10:end))
hold on 
plot(xori(2:10:end))
plot(xori(3:10:end))
plot(xori(4:10:end))
plot(xori(5:10:end))
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

QQ = Baug'*Qaug*Baug+R.*0;
ff = ((Aaug*xR(:,1)-xori)'*Qaug*Baug)';
self.QQ = QQ;
self.ff = ff;



self.sys_info.alpha = 1/max(svd(QQ));


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
% Solver setup
sys_info.alpha = 1/max(svd(QQ)); 
sys_info.epsilon_O = 1e-6;
self.sys_info = sys_info;
self.obs = obs_;
self.dist_arm_id = @(th,DH,base,obs,cap,id) dist_link_Heu(th,DH,base,obs,cap,id);
       
ROBOT = 'M16iB';  
self = CHOMP_FANUC(obs_,sys_info,uref,ROBOT);
self = self.optimizer();
%% solver
% self.nn = horizon*nu;
% self.iter_O = 1;
% self.u = uref;
% self.x_ = x_;
% self.MAX_O_ITER = 10;
% self.cost_all= [];
% tic
% while stop_outer(self)~=true
%    %reset
%    self.cost_old = 1000000;
%    self.cost_new = get_cost_arm(self)+fobs_m(self);%fobs_m(self.sys_info.robot.DH,self.sys_info.robot.base,self.obs{2}.l,self.sys_info.robot,self.x_,self.sys_info.H,self.sys_info.nstate,self.sys_info.njoint,self.obs{2}.D,self.obs{2}.epsilon);
%    self.cost_all = [self.cost_all self.cost_all];
%    
%    %get obstacle cost gradient
%    self.dcobsx = dcostObs_f(self);%dfobs_m(self.sys_info.robot.DH,self.sys_info.robot.base,self.obs{2}.l,self.sys_info.robot,self.x_,self.sys_info.H,self.sys_info.nstate,self.sys_info.njoint,self.obs{2}.D,self.obs{2}.epsilon,self.sys_info.Baug);    
%   
%    %do CHOMP
%    self.cost_old = self.cost_new;
%    self = CHOMP_update_arm(self);
%    %new ref
%    self.x_=[];
%    for i=2:self.sys_info.H+1
%        self.sys_info.xR(:,i)=self.sys_info.robot.A([1:self.sys_info.njoint,7:6+self.sys_info.njoint],[1:self.sys_info.njoint,7:6+self.sys_info.njoint])*self.sys_info.xR(:,i-1)...
%            +self.sys_info.robot.B([1:self.sys_info.njoint,7:6+self.sys_info.njoint],1:self.sys_info.nu)*self.u((i-2)*self.sys_info.nu+1:(i-1)*self.sys_info.nu);
%        self.x_=[self.x_;self.sys_info.xR(:,i)];
%    end   
%    self.iter_O = self.iter_O+1; 
% end   
% toc

uref = self.u;
x_ = self.x_;
cost_all = self.cost_all;

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
[all_ee] = plot_manipulator(ROBOT,obs_{2}.l,xref_,robot,njoint,obs_{2}.D,horizon,robot.base');

%% plot cost all
figure(4)
plot(cost_all)
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Cost')



%% functions

function self = CHOMP_update_arm(self)
   u_ = self.u;
   self.u = u_ - self.sys_info.alpha*3*(dcostArm_f(self) + 100*self.dcobsx);   
  
end

function df = dcostArm_f(self)
           df = self.sys_info.QQ*self.u + self.sys_info.ff;
end

function cost = get_cost_arm(self)
   cost = self.u'*self.QQ*self.u + self.ff'*self.u;
end 
      
function [STOP_O] = stop_outer(self)
   STOP_O  = false;           
   if self.iter_O > self.MAX_O_ITER
       STOP_O  = true;
   end                    
end


function c_all = fobs_m(self)
           DH = self.sys_info.robot.DH;
           base = self.sys_info.robot.base;
           nstate = self.sys_info.nstate;
           njoint = self.sys_info.njoint;  
           c_all = 0;
           for i=1:self.sys_info.H
               theta=self.x_(nstate*(i-1)+1:nstate*(i-1)+njoint);    
               for j = 1:self.obs{1}.num_obs
               Dfx = dm_f(self,theta,j);%dm_f(theta,DH(1:njoint,:),base,self.obs{j+1}.l,self.sys_info.robot.cap,self.obs{j+1}.D);             
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

% 
function dc_all = dcostObs_f(self)
   DH = self.sys_info.robot.DH;
   base = self.sys_info.robot.base;
   nstate = self.sys_info.nstate;
   njoint = self.sys_info.njoint;
   dc_all = zeros(self.nn,1);
   for i=1:self.sys_info.H
       theta=self.x_(nstate*(i-1)+1:nstate*(i-1)+njoint);    
       for j = 1:self.obs{1}.num_obs
           %Dfx = self.dm_f(theta,j);
           Dfx = dm_f(self,theta,j);%dm_f(theta,DH(1:njoint,:),base,self.obs{j+1}.l,self.sys_info.robot.cap,self.obs{j+1}.D);
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





