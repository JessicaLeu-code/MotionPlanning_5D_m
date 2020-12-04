%% main PGDCFS
clear all 
close all
%% build cost

robot=robotproperty(3);
njoint=5;nstate=10;nu=5;DH=robot.DH;
base=[3250, 8500,0]./1000;

load('data/M16_ref_2.mat')
load('data/good_xori.mat')
xref = xuori;
xR=[];xR(:,1)=xref(1:nstate);
horizon=size(uref,1)/njoint;
xori=xref(nstate+1:end);
xref=xori;
uori = uref;
%obs=[[4506;8513;1072]./1000 [4506;8513;1538]./1000];%
obs=[[4106;8313;1]./1000 [4106;8313;1338]./1000];%
%obs=[[4106;8313;1472]./1000 [4106;8313;3038]./1000];%
%% plot ref
xori_ = reshape(xori,[10,24]);
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
D=0.2;
self.epsilon_I = 1e-3;
self.epsilon_O = 1e-4;
self.sys_info.alpha = 1/max(svd(QQ));
%% solver
self.nn = horizon*nu;
self.iter_O = 1;
self.u = uref;
self.MAX_O_ITER = 100;
self.MAX_I_ITER = 20;
self.options = optimoptions('quadprog','Display','off');
self.total_iter = 0;
self.cost_old = 1000000;
self.cost_new = get_cost_arm(self);
cost_all = [];
tic
while stop_outer(self)~=true
   % reset
   self.cost_old = 1000000;
   self.cost_new = get_cost_arm(self);
   self.iter_I = 1;
   
   % get FS
   [self.Ainq,self.binq,Link] = get_con(DH,base,obs,robot,xref,uref,horizon,nstate,njoint,D,nu,Baug);
   % do PSG
   self = inner_PSG_5(self);
   % new ref
   xref=[];
    for i=2:horizon+1
        xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*self.u((i-2)*nu+1:(i-1)*nu);
        xref=[xref;xR(:,i)];
    end
    uref = self.u;
    
    self.QQ = QQ;
    self.ff = ff;    
    cost1 = self.u'*self.QQ*self.u + self.ff'*self.u;
    cost_all = [cost_all cost1];
   self.iter_O = self.iter_O+1; 
end   
toc
%% Plot result

unew = self.u;
xref=[];
for i=2:horizon+1
    xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*unew((i-2)*nu+1:(i-1)*nu);
    xref=[xref;xR(:,i)];
end

xref = reshape(xref,[10,24]);
figure(3)
plot(xref(1:5,:)')
legend('1','2','3','4','5')

uref_ = reshape(self.u,[5,24]);
uori_ = reshape(uori,[5,24]);
figure(2);
plot(uref_','--')
legend('1','2','3','4','5')
hold on;
plot(uori_')

figure(4)
plot(cost_all)


obs = obs.*1000;
%% Visualization (This may take a while ...)

currfold = pwd;
addpath(strcat(currfold,'/Lib'));
addpath(strcat(currfold,'/DERIVESTsuite'));

% xref = xori;
DrawMap;
pause
r = D*700;
[X,Y,Z] = cylinder(r);
hold on
surf(X+obs(1,1),Y+obs(2,1),Z*(obs(3,2)-obs(3,1))+obs(3,1));
[X,Y,Z] = sphere;
surf(X*r+obs(1,1),Y*r+obs(2,1),Z*r+obs(3,1));
surf(X*r+obs(1,1),Y*r+obs(2,1),Z*r+obs(3,2));

plot3(obs(1,1),obs(2,1),obs(3,1),'*')
plot3(obs(1,2),obs(2,2),obs(3,2),'*')
pause

for i=1:1:horizon
    pause(0.1)
    delete(h);
    robot.DH(1:njoint,1)=xref((i-1)*nstate+1:(i-1)*nstate+njoint);
        
    % Use the following lines to draw capsules
%     color=[i/horizon,i/horizon,i/horizon];
%     valpha=0.2;
%     RobotCapLink;
    
    % Use the following lines to draw robot arm
    valpha=1;%i/horizon;
    RobotFigureLink;
    
end



%% functions

function self = inner_PSG_5(self)           
   while stop_inner(self)~=true
       self.cost_old = self.cost_new;
       self = PSG_update_arm(self);
   end
   

end
function self = PSG_update_arm(self)
   u_ = self.u;
   u_ = u_ - self.sys_info.alpha*(dcostArm_f(u_,self.QQ,self.ff) + 1*normrnd(0,0.1,[self.nn,1])./((self.iter_I)^2 + 1));
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
   if delta < self.epsilon_I || self.iter_I > self.MAX_I_ITER
       STOP_I  = true;
   end                    
end
       
function [STOP_O] = stop_outer(self)
   STOP_O  = false;           
   if self.iter_O > self.MAX_O_ITER || norm(self.cost_old - self.cost_new)<self.epsilon_O
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
            %distance=dist_arm_3D(theta,DH(1:njoint,:),base,obs,robot.cap);
            Link = [Link linkid];
            I = [I;distance-D];
            Diff=zeros(njoint,1);
            for s=1:njoint
                [Diff(s),~]=derivest(@(x) dist_link_Heu([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap,linkid),theta(s),'Vectorized','no');
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

