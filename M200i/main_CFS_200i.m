%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFS on manipulator 
% Robot model: FANUC M200i
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

currfold = pwd;
addpath(strcat(currfold,'/Lib'));
addpath(strcat(currfold,'/DERIVESTsuite'));

%% Load robot info
ROBOT = 'M200i';
robot=robotproperty2(ROBOT);
njoint=5;                          % first five joint angles
nstate=10;                         % [ th x 5 ; w x 5 ]
nu=5;                              % [ alpha x 5 ]
DH=robot.DH;


%% Initialization   
load('data/200i_xori.mat')          
x0 = route_wp(:,1);                 % get starting th's 
xg = route_wp(:,end);               % get end th's
xR=[];
xR(:,1)=[x0; zeros(5,1)];
horizon=25;                         % specify horizon/ time duration
% Generatie reference
% line in the state space
xref = [];
for i=1:nstate/2
    xref = [xref; linspace(x0(i),xg(i),horizon+1)];
end%
dt = 0.05; 
for i=1:nstate/2
    w = (xg(i)-x0(i))/(dt*horizon);
    xref = [xref; [w*ones(1,horizon) 0 ]];
end%
xref_ = [];
for i = 1:horizon+1
    xref_ = [xref_; xref(:,i)];
end%  
xori=xref_(nstate+1:end);
xref = xori;
% zeros input
uref = zeros(horizon*njoint,1);
uori = zeros(horizon*njoint,1);

%% Environment
% Set up obstacles
obs=[[3606;8313;1]./1000 [3606;8313;738]./1000];%
D=0.2;                                % minimum distance                                   

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

%% CFS Iteration
self.MAX_O_ITER = 20;
self.iter_O = 0;
total_iter = 0;
tic
cost_all = [];
while stop_outer(self)~=true
    
    % get convex feasible set
    [Lstack,Sstack] = get_con(DH,robot.base,obs,robot,xref,uref,horizon,nstate,njoint,D,nu,Baug);
    
    % solve QP
    options = optimoptions('quadprog','Display','off');
    [unew ,cost,exitflag,output,lambda]= quadprog(Baug'*Qaug*Baug+R.*0,((Aaug*xR(:,1)-xori)'*Qaug*Baug)',Lstack,Sstack,[],[],-0.3*ones(horizon*nu,1),0.3*ones(horizon*nu,1),uref,options);
    self.u = unew;
    uref=unew;
    
    % get cost
    cost1 = self.u'*self.QQ*self.u + self.ff'*self.u;
    cost_all = [cost_all cost1];

    % get new xref
    oldref=xref;
    xref=[];
    for i=2:horizon+1
        xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*unew((i-2)*nu+1:(i-1)*nu);
        xref=[xref;xR(:,i)];
    end    
    norm(xref-oldref)
    % check for termination
    if norm(xref-oldref)<0.01
        disp(strcat('Converged at step',num2str(self.iter_O)));
        break;
    end
    if self.iter_O==self.MAX_O_ITER
         disp('MAX_ITER');
        break;
    end
    
    % next step
    total_iter = total_iter + output.iterations;
    self.iter_O = self.iter_O+1;

end
toc

%% Visualization 
% plot inputs
uref_ = reshape(uref,[5,horizon]);
uori_ = reshape(uori,[5,horizon]);
figure(2);
plot(uref_','--')
legend('\alpha_1','\alpha_2','\alpha_3','\alpha_4','\alpha_5')
xlabel(['Time steps ' '[' num2str(dt) 's/step]'])
ylabel('Input [rad/s^2]')

% plot positions
xref_ = reshape(xref,[10,horizon]);
figure(3);
plot(xref_(1:5,:)','--')
hold on
% plot angle references
plot(xori(1:10:end))
hold on 
plot(xori(2:10:end))
plot(xori(3:10:end))
plot(xori(4:10:end))
plot(xori(5:10:end))
legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5', ...
        '\theta_{1,ref}','\theta_{2,ref}','\theta_{3,ref}','\theta_{4,ref}','\theta_{5,ref}')
xlabel(['Time steps ' '[' num2str(dt) 's/step]'])
ylabel('Input [rad]')

% robot animation
xref_ = reshape(xref,[10,horizon]);
figure(1)
DrawMap2;
delete(h)
[all_ee] = plot_FANUC_200i(obs,xref_,robot,njoint,D,horizon,robot.base');

%% plot cost all
figure(4)
plot(cost_all)
xlabel(['Time steps ' '[' num2str(dt) 's/step]'])
ylabel('Cost')

%% functions

function [Lstack,Sstack] = get_con(DH,base,obs,robot,xref,uref,horizon,nstate,njoint,D,nu,Baug)
    Lstack=[];Sstack=[];
    % for j = 1:num_obs
        I=[];
        for i=1:horizon
            theta=xref(nstate*(i-1)+1:nstate*(i-1)+njoint);
            [distance,linkid]=dist_arm_3D_200i(theta,DH(1:njoint,:),base,obs,robot.cap);
            %distance=dist_arm_3D(theta,DH(1:njoint,:),base,obs,robot.cap);
            I = [I;distance-D];
            Diff=zeros(njoint,1);
            for s=1:njoint
                [Diff(s),~]=derivest(@(x) dist_link_200i([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap,linkid),theta(s),'Vectorized','no');
                %[Diff(s),~]=derivest(@(x) dist_arm_3D([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap),theta(s),'Vectorized','no');
            end            
            %Hess=hessian(@(x) dist_link_Heu(x,DH(1:njoint,:),base,obs,robot.cap,linkid),theta);
            Bj=Baug((i-1)*nstate+1:i*nstate,1:horizon*nu);
            s=I(i)-Diff'*Bj(1:njoint,:)*uref;
            l=-Diff'*Bj(1:njoint,:);
            
%             [E,lambda]=eig(Hess);
%             for m=1:size(lambda,1)
%                 if lambda(m,m)<0
%                     s = s+lambda(m,m)/size(lambda,1);
%                     flag = 'yes';
%                 end
%             end
            
            Sstack=[Sstack;s];
            Lstack=[Lstack;l];
            
        end
    %end
end

%%
function cost = get_cost_arm_f(self)
   cost = self.u'*self.QQ*self.u + self.ff'*self.u;
end 

function [STOP_O] = stop_outer(self)
   STOP_O  = false;           
   if self.iter_O > self.MAX_O_ITER
       STOP_O  = true;
   end                    
end