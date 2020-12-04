%% main CHOMP
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
x0 = xref(1:nstate);
horizon=size(uref,1)/njoint;
xori=xref(nstate+1:end);
xref=xori;
uori = uref;
%obs=[[4506;8513;1072]./1000 [4506;8513;1538]./1000];%
obs=[[4106;8313;1]./1000 [4106;8313;1338]./1000];%
%obs=[[4106;8313;1472]./1000 [4106;8313;3038]./1000];%
D=0.15;
epsilon =  0.05;
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

% constant cost
caug = (Aaug*xR(:,1)-xori)'*Qaug*(Aaug*xR(:,1)-xori);

sys_info.Aaug = Aaug;
sys_info.Baug = Baug;
sys_info.Qaug = QQ;
sys_info.paug = ff;
sys_info.caug = caug;

self.QQ = QQ;
self.ff = ff;


self.epsilon = 1e-6; 
self.sys_info.alpha = 1/max(svd(QQ));
%% solver
self.nn = horizon*nu;
self.iter_O = 1;
self.u = uref;
self.MAX_O_ITER = 10;
sys_info.lambda = 10000;
self.options = optimoptions('quadprog','Display','off');
% tic
% while stop_outer(self)~=true
%    % reset
%    self.cost_old = 1000000;
%    self.cost_new = get_cost_arm(self)+fobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon);
%    self.iter_I = 1;
%    
%    % get obstacle cost gradient
%    self.dcobsx = dfobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon,Baug);    
%   
%    % do CHOMP
%    self.cost_old = self.cost_new;
%    self = CHOMP_update_arm(self);
%    % new ref
%    xref=[];
%     for i=2:horizon+1
%         xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*self.u((i-2)*nu+1:(i-1)*nu);
%         xref=[xref;xR(:,i)];
%     end
%     uref = self.u;
%    self.iter_O = self.iter_O+1; 
% end   
% toc

%HMC
%self.u = zeros(horizon*5,1);
% solver
ee = 0.00025;
u_all = [];
cost_all = [];
pp_all =[];
tic

%% test values





%% solver
% for k = 1:1
%     
%     % resample
%     ut = self.u;
%     %rt = (rand(njoint*horizon,1)-0.5);%*(self.sys_info.alpha/k);
%     rt = normrnd(0,1,[njoint*horizon,1]);
%     
%     for kk = 1:100
%         
%         kk
%         % leapfrog
%         L = 10;
%         old_rt = rt;
%         old_xref =  xref;
%         rt_ = rt - (ee/2)*dU_f(self,sys_info,DH,base,obs,robot,uref,xR,horizon,nstate,njoint,D,epsilon,Baug); 
%         old_u = self.u;
%         for ll = 1:L
%             self.u = ut + ee*rt_;
%             rt_ = rt_ -(ee)*dU_f(self,sys_info,DH,base,obs,robot,self.u,xR,horizon,nstate,njoint,D,epsilon,Baug);
%         end
%         rt_ee = rt_ -(ee/2)*dU_f(self,sys_info,DH,base,obs,robot,self.u,xR,horizon,nstate,njoint,D,epsilon,Baug);
%         rt_ee = -rt_ee;
%         
%         xref=[];
%         for i=2:horizon+1
%             xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*self.u((i-2)*nu+1:(i-1)*nu);
%             xref=[xref;xR(:,i)];
%         end
%         uref = self.u;
%         
%         if H_f(sys_info,obs,rt,old_u,k,DH,base,robot,old_xref,horizon,nstate,njoint,D,epsilon)<H_f(sys_info,obs,rt_ee,self.u,k,DH,base,robot,xref,horizon,nstate,njoint,D,epsilon)
%             rt = rt_ee;
%             ut = self.u;
%             pp_all =[pp_all 1];
%         else
%             pp = H_f(sys_info,obs,rt_ee,self.u,k,DH,base,robot,xref,horizon,nstate,njoint,D,epsilon)/H_f(sys_info,obs,rt,old_u,k,DH,base,robot,old_xref,horizon,nstate,njoint,D,epsilon);
%             pp_all =[ pp_all pp];
%             if rand<pp
%                 rt = rt_ee;
%                 ut = self.u;
%             else
%                 self.u = old_u;
% %                 rt = (rand(njoint*horizon,1)-0.5);%*(self.sys_info.alpha);
% %                 rt = norm(old_rt)*rt/norm(rt);
%                   rt = normrnd(0,1,[njoint*horizon,1]);
%             end
%         end
%         u_all = [ u_all self.u];
%         cost_all = [cost_all getCost(self,sys_info,obs,DH,base,robot,xref,horizon,nstate,njoint,D,epsilon)];
%         % new ref
%         xref=[];
%         for i=2:horizon+1
%             xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*self.u((i-2)*nu+1:(i-1)*nu);
%             xref=[xref;xR(:,i)];
%         end
%         uref = self.u;
%     end
%     
% end
[lpdf,glpdf] = LPDFGLPDF(self.u,self,sys_info,DH,base,obs,robot,x0,horizon,nstate,njoint,D,epsilon,Baug)


logpdf = @(u)LPDFGLPDF(u,self,sys_info,DH,base,obs,robot,x0,horizon,nstate,njoint,D,epsilon,Baug);
smp = hmcSampler(logpdf,self.u);
smp = tuneSampler(smp);
u_all = drawSamples(smp,'NumSamples',100,'VerbosityLevel',1);

toc
%%
cost_all =[];
for i = 1:size(u_all,1)
     cost_all = [cost_all U_f2(u_all(i,:)',sys_info,obs)];     
end


%% Cluster
numc = 10;
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
'Streams',stream);
tic; % Start stopwatch timer
[idx,C,sumd,D] = kmeans(u_all,numc,'Options',options,'MaxIter',10000,...
'Display','final','Replicates',10);toc;

%%
figure
plot(idx)
hold on
plot(cost_all)



%%

figure 
plot(pp_all,'b')

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

function self = CHOMP_update_arm(self)
   u_ = self.u;
   self.u = u_ - self.sys_info.alpha*3*(dcostArm_f(u_,self.QQ,self.ff) + 100*self.dcobsx);   
  
end

function df = dcostArm_f(u_,H,f)
    df = 2*(H*u_ + f);
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

function c_all = fobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon)
c_all = 0;
for i=1:horizon
    theta=xref(nstate*(i-1)+1:nstate*(i-1)+njoint);    
    Dfx = dm_f(theta,DH(1:njoint,:),base,obs,robot.cap,D);        
    for s =1:njoint
        if Dfx(s)<0
            c_x = -Dfx(s)+(1/2)*epsilon;    
        elseif Dfx(s)>=0 && Dfx(s)<=epsilon
            c_x = (1/(2*epsilon))*(Dfx(s)-epsilon)^2;
        else
            c_x = 0;
        end

        c_all = c_all + c_x;
    end
end

end

function dc_all = dfobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon,Baug)
dc_all = zeros(horizon*njoint,1);
for i=1:horizon
    theta=xref(nstate*(i-1)+1:nstate*(i-1)+njoint);    
    Dfx = dm_f(theta,DH(1:njoint,:),base,obs,robot.cap,D);
    [dis, linkid] = min(Dfx);
    dDfx = zeros(njoint,1);
    Bj=Baug((i-1)*nstate+1:i*nstate,:);
    
    if Dfx(linkid)<0
        for s=1:njoint
            [dDfx(s),~]=derivest(@(x) dist_link_Heu([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap,linkid),theta(s),'Vectorized','no');
            %[Diff(s),~]=derivest(@(x) dist_arm_3D([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap),theta(s),'Vectorized','no');
        end            
        dc_x = -[dDfx'*Baug((i-1)*njoint+1:i*njoint,:)]';

    elseif Dfx(linkid)>=0 && Dfx(linkid)<=epsilon
        for s=1:njoint
            [dDfx(s),~]=derivest(@(x) dist_link_Heu([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap,linkid),theta(s),'Vectorized','no');
            %[Diff(s),~]=derivest(@(x) dist_arm_3D([theta(1:s-1);x;theta(s+1:end)],DH(1:njoint,:),base,obs,robot.cap),theta(s),'Vectorized','no');
        end 

        dc_x = [(1/epsilon)*(Dfx(linkid)-epsilon)*dDfx'*Baug((i-1)*njoint+1:i*njoint,:)]';
    else
        dc_x =  zeros(horizon*njoint,1);
    end

    dc_all  = dc_all + dc_x;
    
    
end

end

function  [dU] = dU_f(self,sys_info,DH,base,obs,robot,uref,xR,horizon,nstate,njoint,D,epsilon,Baug)
        xref=[];
        for i=2:horizon+1
            xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:njoint)*self.u((i-2)*njoint+1:(i-1)*njoint);
            xref=[xref;xR(:,i)];
        end
        self.dcobsx = dfobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon,Baug);
        dU = (dcostArm_f(uref,self.QQ,self.ff)+sys_info.lambda*self.dcobsx);
end

function [Hfx] = H_f(sys_info,obs,rt,u,k,DH,base,robot,xref,horizon,nstate,njoint,D,epsilon)
    
    cobs = fobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon);
    Ux = 0.5*u'*sys_info.Qaug*u +sys_info.paug'*u + sys_info.caug + cobs;
    Kx = 0.5*rt'*rt;
    Hfx = exp(-(k/2000)*(Ux+Kx));
end

function Ux = getCost(self,sys_info,obs,DH,base,robot,xref,horizon,nstate,njoint,D,epsilon)
    cobs = fobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon);
    Ux = 0.5*self.u'*sys_info.Qaug*self.u +sys_info.paug'*self.u + sys_info.caug + sys_info.lambda*cobs;
end

% matlab HMC

function xref_ = get_xref(u,horizon,njoint,robot,x0)
    xref_=[];
    xR(:,1) = x0;
    for i=2:horizon+1
        xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:njoint)*u((i-2)*njoint+1:(i-1)*njoint);
        xref_=[xref_;xR(:,i)];
    end
end

function [lpdf,glpdf] = LPDFGLPDF(u,self,sys_info,DH,base,obs,robot,x0,horizon,nstate,njoint,D,epsilon,Baug)
xref = get_xref(u,horizon,njoint,robot,x0);
lpdf = -U_f2(u,sys_info,obs,DH,base,robot,xref,horizon,nstate,njoint,D,epsilon);
glpdf = -dU_f2(u,self,sys_info,DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon,Baug);
end

function Ux = U_f2(u,sys_info,obs,DH,base,robot,xref,horizon,nstate,njoint,D,epsilon)
    cobs = fobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon);
    Ux = 0.5*u'*sys_info.Qaug*u +sys_info.paug'*u + sys_info.caug + sys_info.lambda*cobs;
end

function  [dU] = dU_f2(u,self,sys_info,DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon,Baug)
        
        self.dcobsx = dfobs_m(DH,base,obs,robot,xref,horizon,nstate,njoint,D,epsilon,Baug);
        dU = (dcostArm_f(u,self.QQ,self.ff)+sys_info.lambda*self.dcobsx);
end


