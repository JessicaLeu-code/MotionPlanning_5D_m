%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFS on robot arm 
% In this example, the arm is FANUC M16iB
%
% Changliu Liu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

currfold = pwd;
addpath(strcat(currfold,'/Lib'));
addpath(strcat(currfold,'/DERIVESTsuite'));

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
%obs=[[4106;8313;1]./1000 [4106;8313;2238]./1000];%
%obs=[[4106;8313;1]./1000 [4106;8313;2238]./1000];%
obs=[[4106;8313;1]./1000 [4106;8313;1338]./1000];%
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
%% CFS Iteration
MAX_ITER = 10;
total_iter = 0;
D=0.2;
tic
cost_all = [];
for k=1:MAX_ITER

    [Lstack,Sstack] = get_con(DH,base,obs,robot,xref,uref,horizon,nstate,njoint,D,nu,Baug);

    %pretime=toc;
    %tic;
    options = optimoptions('quadprog','Display','off');
    [unew ,cost,exitflag,output,lambda]= quadprog(Baug'*Qaug*Baug+R.*0,((Aaug*xR(:,1)-xori)'*Qaug*Baug)',Lstack,Sstack,[],[],-0.2*ones(horizon*nu,1),0.2*ones(horizon*nu,1),uref,options);
    total_iter = total_iter + output.iterations;
    QQ = Baug'*Qaug*Baug+R.*0;
    ff = ((Aaug*xR(:,1)-xori)'*Qaug*Baug)';
    self.QQ = QQ;
    self.ff = ff;    
    self.u = unew;
    cost1 = self.u'*self.QQ*self.u + self.ff'*self.u;
    cost_all = [cost_all cost1];
    %     time=toc
%     
%     fun= @(u) u'* (Baug'*Qaug*Baug+R)* u+ 2*((Aaug*xR(:,1)-xori)'*Qaug*Baug)*u;%+(Aaug*xR(:,1)-xref)'*Qaug*(Aaug*xR(:,1)-xref);
%     nonlcon=@(u) nldisfn(u,robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint]),robot.B([1:njoint,7:6+njoint],1:nu),xori(1:nstate),robot.DH,horizon,base,obs,robot.cap);
%     options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','TolCon',1,'MaxFunEvals',30000);
%     tic;
%     [unew,cost2 ]= fmincon(fun,uref,[],[],[],[],[],[],[],options);
%     time=toc
    
    oldref=xref;
    xref=[];
    for i=2:horizon+1
        xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*unew((i-2)*nu+1:(i-1)*nu);
        xref=[xref;xR(:,i)];
    end
    uref=unew;
    
    if norm(xref-oldref)<0.01
        disp(strcat('Converged at step',num2str(k)));
        break;
    end
    if k==MAX_ITER
         disp('MAX_ITER');
        break;
    end
%     if Lstack*uref<Sstack
%         disp(strcat('Local Optimal Found at step ',num2str(k)));
%         break;
%     end
end
toc
%%
% unew = uori;
% xref=[];
%     for i=2:horizon+1
%         xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*unew((i-2)*nu+1:(i-1)*nu);
%         xref=[xref;xR(:,i)];
%     end
%%
% unew = uori;
% xuori=[];
%     for i=2:horizon+1
%         xR(:,i)=robot.A([1:njoint,7:6+njoint],[1:njoint,7:6+njoint])*xR(:,i-1)+robot.B([1:njoint,7:6+njoint],1:nu)*unew((i-2)*nu+1:(i-1)*nu);
%         xuori=[xuori;xR(:,i)];
%     end
% xuori_ = reshape(xuori,[10,24]);
%  %
% plot(xuori_(1:5,:)')
% legend('1','2','3','4','5')

QQ = Baug'*Qaug*Baug+R.*0;
ff = ((Aaug*xR(:,1)-xori)'*Qaug*Baug)';
self.QQ = QQ;
self.ff = ff;
self.u = unew;
cost = self.u'*self.QQ*self.u + self.ff'*self.u;
   
obs = obs.*1000;
%% Visualization (This may take a while ...)
%xref = xori;
DrawMap;
%%
[X,Y,Z] = cylinder(D*800);
hold on
surf(X+obs(1,1),Y+obs(2,1),Z*(obs(3,2)-obs(3,1))+obs(3,1));
[X,Y,Z] = sphere;
surf(X*D*800+obs(1,1),Y*D*800+obs(2,1),Z*D*800+obs(3,1));
surf(X*D*800+obs(1,1),Y*D*800+obs(2,1),Z*D*800+obs(3,2));
plot3(obs(1,1),obs(2,1),obs(3,1),'*')
plot3(obs(1,2),obs(2,2),obs(3,2),'*')
%%
for i=1:1:horizon
    robot.DH(1:njoint,1)=xref((i-1)*nstate+1:(i-1)*nstate+njoint);
        
    % Use the following lines to draw capsules
%     color=[i/horizon,i/horizon,i/horizon];
%     valpha=0.2;
%     RobotCapLink;
    
    % Use the following lines to draw robot arm
    valpha=1;%i/horizon;
    RobotFigureLink;

end
%%
uref_ = reshape(uref,[5,24]);
uori_ = reshape(uori,[5,24]);
figure(2);
plot(uref_','--')
legend('1','2','3','4','5')
hold on;
plot(uori_')

xref_ = reshape(xref,[10,24]);
figure(3);
plot(xref_(1:5,:)','--')
legend('1','2','3','4','5')

%% plot cost all
figure(4)
plot(cost_all)


function [Lstack,Sstack] = get_con(DH,base,obs,robot,xref,uref,horizon,nstate,njoint,D,nu,Baug)
    Lstack=[];Sstack=[];
    % for j = 1:num_obs
        I=[];
        for i=1:horizon
            theta=xref(nstate*(i-1)+1:nstate*(i-1)+njoint);
            [distance,linkid]=dist_arm_3D_Heu(theta,DH(1:njoint,:),base,obs,robot.cap);
            %distance=dist_arm_3D(theta,DH(1:njoint,:),base,obs,robot.cap);
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

%%
function cost = get_cost_arm_f(self)
   cost = self.u'*self.QQ*self.u + self.ff'*self.u;
end 