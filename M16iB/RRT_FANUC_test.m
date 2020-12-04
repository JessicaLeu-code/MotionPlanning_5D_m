%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RRT and RRT* on robot arm 
% In this example, the arm is FANUC M16iB
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

%% Environment
% Obstacle
%obs=[[4106;8313;1]./1000 [4106;8313;1938]./1000];
obs=[[3806;8313;1]./1000 [3806;8313;1938]./1000];
num_obs = 1;
obs_info.num_obs = num_obs;
obs_info.obs = obs;
obs_info.D=0.2;
offset = [3250 8500 0];

horizon = 1;
%% Robot
currfold = pwd;
addpath(strcat(currfold,'/Lib'));
addpath(strcat(currfold,'/DERIVESTsuite'));

robot=robotproperty(3);
njoint=5;nstate=10;nu=5;DH=robot.DH;
base=[3150, 8500,0]./1000;

load('data/M16_ref_2.mat')
load('data/good_xori.mat')

% initialization 
% number of states
nstate = 5;
x0 = [0.4921, 1.5922, -0.0092, -0.0010, -1.5786]';

    
% Goal
goalxyz = [-0.6090 1.0873 0.3008 -0.0515 -1.4130]';
goal_th = [-0.6090 1.0873 0.3008 -0.0515 -1.4130]';
region_g = [pi/20   pi/20  pi/10 pi/2  pi/2]';


%% State constraints
region_s = [pi/2   pi/2  pi/2    pi/1.5  pi/1.5]';
sample_off = [0 pi/2  0 0 0]';
ratial = [1 1 0.5 0.1 0.1]';

%% System info
sys_info.DH = DH;
sys_info.robot = robot;
sys_info.nstate = nstate;
sys_info.x0 = x0;
sys_info.base = base;
sys_info.ratial = ratial;
% plot initial 
figure(1)
DrawMap;
delete(h)
plot_FANUC(obs,x0,robot,njoint,obs_info.D,1,offset);
sys_info.goal_th = goal_th;
sys_info.fighandle = figure(1);   


%% solver RRT

tic
self = RRT_FANUC(obs_info,sys_info,goalxyz,region_g,region_s,sample_off);
%self = self.RRTfind();
self = self.RRTstar_find();
Time_{2}= toc 

%% solver spmd RRT

% 
% pool_num = 5;
% parpool(pool_num);
% tic
% spmd
%     
%     self = RRT_FANUC(obs_info,sys_info,goalxyz,region_g,region_s,sample_off);
%     % initialize tree
%     self.newNode = [self.sys_info.x0];
%     self.all_nodes = [-1; self.newNode];
%     self = self.goal_reached();
%     
%     done = false;
%     while ~done
%         self = self.getNode();
%         self = self.addNode();
%         self = self.goal_reached();
%         done = gop(@or, self.isGoalReached);
%     end 
%            
%     self.route = self.newNode;
%     while self.parent ~= -1
%         self.route = [ self.all_nodes(2:end,self.parent) self.route];
%         self.parent = self.all_nodes(1,self.parent);
%     end
%     
% end
% 
% for j = 1:pool_num
%     cc = self{j};
%     if cc.isGoalReached == 1
%         route_id = j;
%         break;
%     end
% end
% 
% cc.route = cc.newNode;
% cc.parent = cc.all_nodes(1,end);
% while cc.parent ~= -1
%     cc.route = [ cc.all_nodes(2:end,cc.parent) cc.route];
%     cc.parent = cc.all_nodes(1,cc.parent);
% end
% 
% 
% Time_{2}= toc 
% 
% self = cc;
%% Plot
path_length = size(self.route,2);
figure(1)
DrawMap;
delete(h)
[all_ee] = plot_FANUC(obs,self.route,robot,njoint,obs_info.D,path_length,offset);


%delete(gcp('nocreate'))

cc = self;
cc.all_ee = cc.all_ee.*1000;
plot3(cc.all_ee(1,:),cc.all_ee(2,:),cc.all_ee(3,:),'g*')


plot3(all_ee(1,:),all_ee(2,:),all_ee(3,:),'ro')


%% sample
num_of_samples = 20;


gap = 3;
route_wp = [];
ee_wp = [];
for  i=1:gap:path_length
    route_wp = [route_wp self.route(:,i)];
    ee_wp = [ee_wp all_ee(:,i)];
end

plot3(ee_wp(1,:),ee_wp(2,:),ee_wp(3,:),'b*')
    



