%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RRT and RRT* on manipulator 
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
goalxyz = [-0.8090 0.1873 0.3008 0.1515 -1.4130]';
goal_th = goalxyz;
region_g = [pi/20 pi/20 pi/10 pi/2 pi/2]';   % goal region

%% Environment
% Set up obstacles
obs=[[3606;8413;1]./1000 [3606;8413;1038]./1000]; % cylinder (straight line)
num_obs = 1;
obs_info.num_obs = num_obs;
obs_info.obs = obs;
obs_info.D=0.2;

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
sys_info.fighandle = figure(1);   


%% solver RRT
tic
self = RRT_FANUC(obs_info,sys_info,goalxyz,region_g,region_s,sample_off,'M200i','RRT*');
self = self.find_route();
Time_{2}= toc 

%% Plot
path_length = size(self.route,2);
% robot animation
figure(2)
DrawMap2;
delete(h)
[all_ee] = plot_FANUC_200i(obs,self.route,robot,njoint,obs_info.D,path_length,offset);

cc = self;
% plot all samples
plot3(cc.all_ee(1,:),cc.all_ee(2,:),cc.all_ee(3,:),'g*')
% plot eee route
plot3(all_ee(1,:),all_ee(2,:),all_ee(3,:),'ro')

%% Sample for CFS 
gap = 2;
route_wp = [];
ee_wp = [];
for  i=1:gap:path_length
    route_wp = [route_wp self.route(:,i)];
    ee_wp = [ee_wp all_ee(:,i)];
end

% plot samples
plot3(ee_wp(1,:),ee_wp(2,:),ee_wp(3,:),'b*')
    



