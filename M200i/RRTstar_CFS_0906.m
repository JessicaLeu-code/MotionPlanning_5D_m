%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RRT*-CFS on manipulator 
% Robot model: FANUC M200i
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

currfold = pwd;
addpath(strcat(currfold,'/Lib'));
addpath(strcat(currfold,'/DERIVESTsuite'));

%% Load robot info
ROBOT = 'M200i';
%load(['figure/',ROBOT,'.mat']);
robot=robotproperty2(ROBOT);
njoint=5;                          % first five joint angles
nstate=10;                         % [ th x 5 ; w x 5 ]
nu=5;                              % [ alpha x 5 ]
DH=robot.DH;
offset= [3150, 8500,330]./1000;    % center of the scene
robot.base = offset';

%% Initialization 
% number of states
nstate = 5;
% initial position
x0 = [0.421, 0, -0.0092, -0.0010, -1.5786]'; % first five angles

%% Goal
goalxyz = [-1.4090 0.8873 0.4008 0.0 0.4430]';
%goalxyz = [-0.7825; 0.0284; 0.2172; 0.1444; -1.1779];  
goal_th = goalxyz;
region_g = [pi/20 pi/20 pi/10 pi/2 pi/2]';   % goal region

%% Environment
% Set up obstacles
obs{1}.shape = 'cylinder';
obs{1}.A = [1 0 ; 0 1];
obs{1}.l = [[3606;8413;1]./1000 [3606;8413;1038]./1000]; % cylinder (straight line);% % for M200i
obs{1}.D = 0.2;
obs{1}.epsilon =  0.2;

obs{2}.shape = 'cylinder';
obs{2}.A = [1 0 ; 0 1];
obs{2}.l = [[3406;7813;800]./1000 [3406;7813;1538]./1000]; % cylinder (straight line);% % for M200i
obs{2}.D = 0.2;
obs{2}.epsilon =  0.2;

%% State constraints
region_s = [pi/2   pi/2  pi/2    pi/1.5  pi/1.5]';
sample_off = [0 0  0 0 0]';
ratial = [1 1 0.5 0.1 0.1]';

%% System info
sys_info.DH = DH;
sys_info.robot = robot;
sys_info.nstate = nstate;
sys_info.x0 = x0;
sys_info.base = robot.base;
sys_info.ratial = ratial;
sys_info.goal_th = goal_th;

%% solver RRT
% % single rrt
% tic
% self = RRT_FANUC(obs,sys_info,goalxyz,region_g,region_s,sample_off,'M200i','RRT*');
% self = self.find_route();
% Time_(1)= toc 

% % multi-thread rrt
tic
%%%%%%%%%%%%%%%%%%%%
s_Parallel_rrt;
%%%%%%%%%%%%%%%%%%%%
Time_(1)= toc; 
%% Plot
path_length = size(self.route,2);
% robot animation
figure(1)
DrawMap2;
delete(h)
[all_ee] = plot_200i(obs,self.route,robot,njoint);

cc = self;
% plot all samples
p1 = [];
p1(1) = plot3(cc.all_ee(1,:),cc.all_ee(2,:),cc.all_ee(3,:),'g*');
% plot eee route
p1(2) = plot3(all_ee(1,:),all_ee(2,:),all_ee(3,:),'co');

%% Sample for CFS 
nstate = 10;
dt = robot.delta_t;
horizon = 40;
wpTimes = (0:size(self.route,2)-1)*dt;
trajTimes = linspace(0,wpTimes(end),horizon+1);
sampled_route = cubicpolytraj(self.route,wpTimes,trajTimes);
figure(1)
[sampled_ee] = plot_200i(obs,sampled_route,robot,njoint);
p1(3) = plot3(sampled_ee(1,:),sampled_ee(2,:),sampled_ee(3,:),'bo-','LineWidth',1);

%% Generate ref for CFS
tic
% Get reference
xref = sampled_route;
x0 = xref(:,1);
xg = xref(:,end);
% add w refernece
xref = [xref;zeros(5,horizon+1)];   % penalize angular velocity
xref_ = [];
for i = 1:horizon+1
    xref_ = [xref_; xref(:,i)];
end
xori=xref_(nstate+1:end);
xR(:,1)=xref_(1:nstate);
xref=xori;
% zeros input  
uref = zeros(horizon*njoint,1);
uori = zeros(horizon*njoint,1);
%% Cost Fn Parameters
Aaug=[];Baug=zeros(horizon*nstate,horizon*nu);Qaug=zeros(horizon*nstate);
Q=[];
Q(1:njoint,1:njoint)=[10 0 0 0 0;
    0 10 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];
Q(1:njoint,njoint+1:2*njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,1:njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,njoint+1:2*njoint)=[100 0 0 0 0;
    0 20 0 0 0;
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
    R((i-1)*nu+1:i*nu,(i-1)*nu+1:i*nu)=[10 0 0 0 0;
        0 10 1 0 0;
        0 1 2 0 0;
        0 0 0 2 0;
        0 0 0 0 1];
end
R=R+R';
QQ = Baug'*Qaug*Baug+R.*10;
gaug = kron(ones(horizon,1),[xg;0;0;0;0;0]);  %don't use. Bad performance.
% not quite.
%gaug = xref;
ff = ((Aaug*xR(:,1)-gaug)'*Qaug*Baug)';
% constant cost
caug = (Aaug*xR(:,1)-gaug)'*Qaug*(Aaug*xR(:,1)-gaug);

%% System info
sys_info.Aaug = Aaug;
sys_info.Baug = Baug;
sys_info.QQ = QQ;
sys_info.ff = ff;
sys_info.Qaug = QQ;
sys_info.paug = ff;
sys_info.caug = caug;
sys_info.robot = robot;
sys_info.H = horizon;
sys_info.nstate = nstate;
sys_info.njoint = njoint;
sys_info.xR = xR;
sys_info.nu = nu;
sys_info.x_ = xori;
sys_info.alpha = 1/max(svd(QQ)); 
sys_info.lim = [1;1;1;1;1];
% Solver setup
sys_info.options = optimoptions('quadprog','Display','off');
sys_info.alpha = 1/max(svd(QQ)); 
sys_info.epsilon_O = 1e-1;
sys_info.MAX_O_ITER = 20;  % 6 for chomp
sys_info.MAX_input = kron(ones(horizon,1),[1;1;pi;pi;pi]*dt);

%% Solver
% baseline
eval = EVAL(sys_info);
Cost_b = eval.get_Cost_b();
% CFS 
self_cfs = CFS_FANUC(obs,sys_info,ROBOT);
self_cfs = self_cfs.optimizer();
Time_(2)= toc
cost_all = self_cfs.eval.cost_new
iter_rrt = 1;
iter_rrt_cfs = self_cfs.iter_O-1+iter_rrt

uref = self_cfs.u;
x_ = self_cfs.x_;
xref_ = reshape(x_,[nstate,horizon]);

%% Visualization 
% robot animation
figure(1)
[all_ee] = plot_manipulator2(ROBOT,obs,xref_,robot,njoint);
axis([2.5 4.5 7 9.5 -0.01 2])
p1(4) = plot3(all_ee(1,:),all_ee(2,:),all_ee(3,:),'ko');
legend(p1,'RRT* samples','RRT* solution','sampled route','CFS solution')

% plot inputs
uref_ = reshape(uref,[5,horizon]);
uori_ = reshape(uori,[5,horizon]);
figure(3);
plot(uref_','--')
legend('\alpha_1','\alpha_2','\alpha_3','\alpha_4','\alpha_5')
xlabel(['Time steps ' '[' num2str(dt) 's/step]'])
ylabel('Input [rad/s^2]')

% plot positions
xref_ = reshape(x_,[10,horizon]);
figure(4);
plot(xref_(1:5,:)')
hold on
% plot angle references
plot(xori(1:10:end),'--')
hold on 
plot(xori(2:10:end),'--')
plot(xori(3:10:end),'--')
plot(xori(4:10:end),'--')
plot(xori(5:10:end),'--')
legend('\theta_1','\theta_2','\theta_3','\theta_4','\theta_5', ...
        '\theta_{1,ref}','\theta_{2,ref}','\theta_{3,ref}','\theta_{4,ref}','\theta_{5,ref}','Location','northoutside','Orientation','horizontal','NumColumns',5);
xlabel(['Time steps ' '[' num2str(dt) 's/step]'])
ylabel('Input [rad]')

%% plot cost all
figure(5)
plot(self_cfs.eval.cost_all)
xlabel(['Time steps ' '[' num2str(dt) 's/step]'])
ylabel('Cost')
