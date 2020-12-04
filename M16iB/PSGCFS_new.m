%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PGDCFS on manipulator 
% Robot model: FANUC M16iB
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

%% Load robot info
ROBOT = 'M16iB';
robot = robotproperty2(ROBOT);
njoint=5;nstate=10;nu=5;DH=robot.DH;
offset = [3250, 8500,0]./1000;   % center of the scene
robot.base = offset';

%% Initialization  
load('data/good_xori.mat')
xR=[];
xR(:,1) = xuori(1:nstate);
x0 = xuori(1:nstate);            % get starting th's 
xg = xuori(end-nstate+1:end);    % get end th's
horizon = 24;                    % specify horizon/ time duration
% Generatie reference
% line in the state space
xref = [];
for i=1:nstate/2
    xref = [xref; linspace(x0(i),xg(i),horizon+1)];
end%
W = zeros(5,horizon+1);
xref = [xref;W];

xref_ = [];
for i = 1:horizon+1%
    xref_ = [xref_; xref(:,i)];
end   
xori=xref_(nstate+1:end);
xref = xori;
% zeros input
uref = zeros(horizon*njoint,1);
uori = zeros(horizon*njoint,1);

%% Environment
% Set up obstacles
obs=[[4106;8313;1]./1000 [4106;8313;1338]./1000];%
D=0.2;

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
QQ = Baug'*Qaug*Baug+R.*10;
ff = ((Aaug*xR(:,1)-xori)'*Qaug*Baug)';
self.QQ = QQ;
self.ff = ff;

%% solver
% Parameter
self.nn = horizon*nu;
% Solver setup
self.sys_info.alpha = 1/max(svd(QQ));    % stochastic ratio 
self.MAX_O_ITER = 10;
self.MAX_I_ITER = 20;
self.epsilon = 1e-6; 
self.options = optimoptions('quadprog','Display','off');
% Initialize
self.iter_O = 0;
self.total_iter = 0;
self.u = uref;
self.u_old = ones(self.nn,1);
self.cost_old = 1000000;
self.cost_new = get_cost_arm(self);
cost_all = [];

tic
while stop_outer(self)~=true
    % reset
    self.u_old = self.u;    
    self.cost_new = get_cost_arm(self);
    self.iter_I = 1;
   
    % get feasible sert
    [self.Ainq,self.binq,Link] = get_con(DH,robot.base,obs,robot,xref,uref,horizon,nstate,njoint,D,nu,Baug);
    
    % do PSG
    self = inner_PSG_5(self);
    
    % get cost
    cost1 = self.u'*self.QQ*self.u + self.ff'*self.u;
    cost_all = [cost_all cost1];
    
    % get new xref
    xref=[];
    for i=2:horizon+1
        xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*self.u((i-2)*nu+1:(i-1)*nu);
        xref=[xref;xR(:,i)];
    end    
    uref = self.u;   
    
    % next step
    self.iter_O = self.iter_O+1;
end 
toc

%% Visualization 
uref_ = reshape(self.u,[5,24]);
uori_ = reshape(uori,[5,24]);
figure(2);
plot(uref_','--')
legend('1','2','3','4','5')
hold on;
plot(uori_')
%% Plot result
% plot inputs
uref_ = reshape(uref,[5,horizon]);
uori_ = reshape(uori,[5,horizon]);
figure(2);
plot(uref_','--')
legend('\alpha_1','\alpha_2','\alpha_3','\alpha_4','\alpha_5')
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Input [rad/s^2]')

% plot positions
xref = reshape(xref,[10,24]);
figure(3)
plot(xref(1:5,:)')
hold on
legend('1','2','3','4','5')
% plot ref
xori_ = reshape(xori,[10,24]);
plot(xori_(1:5,:)','--')
legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5', ...
        '\theta_{1,ref}','\theta_{2,ref}','\theta_{3,ref}','\theta_{4,ref}','\theta_{5,ref}')
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Input [rad]')

% robot animation
xref_ = reshape(xref,[10,horizon]);
figure(1)
DrawMap;
delete(h)
[all_ee] = plot_FANUC(obs,xref_,robot,njoint,D,horizon,offset);

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
   u_ = u_ - self.sys_info.alpha*(dcostArm_f(u_,self.QQ,self.ff) + 0*normrnd(0,0.1,[self.nn,1])./((self.iter_I)^2 + 1));
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

function df = dcostArm_f(u_,H,f)
    df = H*u_ + f;
end

function cost = get_cost_arm(self)
   cost = self.u'*self.QQ*self.u + self.ff'*self.u;
end 

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

function [Lstack,Sstack,Link] = get_con(DH,base,obs,robot,xref,uref,horizon,nstate,njoint,D,nu,Baug)
    Lstack=[];Sstack=[];
    % for j = 1:num_obs
        I=[];
        Link = [];
        for i=1:horizon
            theta=xref(nstate*(i-1)+1:nstate*(i-1)+njoint);
            [distance,linkid]=dist_arm_3D_Heu(theta,DH(1:njoint,:),base,obs,robot.cap);
            Link = [Link linkid];
            I = [I;distance-D];
            Diff=zeros(njoint,1);
            for s=1:njoint
                [Diff(s),~]=derivest(@(x) dist_link_Heu([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap,linkid),theta(s),'Vectorized','no');
            end            
            %Hess=hessian(@(x) dist_link_Heu(x,DH(1:njoint,:),base,obs,robot.cap,linkid),theta);
            Bj=Baug((i-1)*nstate+1:i*nstate,1:horizon*nu);
            s=I(i)-Diff'*Bj(1:njoint,:)*uref;
            l=-Diff'*Bj(1:njoint,:);           
            
            Sstack=[Sstack;s];
            Lstack=[Lstack;l];
            
        end
    %end
end

