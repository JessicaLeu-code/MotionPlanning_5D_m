%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planning comparison on manipulator
% With planning algorithm: CFS, PSGCFS
% Robot model: FANUC M200i, M16iB
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

%% Load robot info
ROBOT = 'M200i';                   % 'M200i' or 'M16iB'
robot = robotproperty2(ROBOT);

njoint=5;                          % first five joint angles
nstate=10;                         % [ th x 5 ; w x 5 ]
nu=5;                              % [ alpha x 5 ]
DH=robot.DH;

%% Solver set up
% SOLVER{2}.name = 'PSGCFS';
% SOLVER{2}.plot = false;
SOLVER{3}.name ='CFS';
SOLVER{3}.plot = true;


%% Initialization  
% % M200i
x0 = [0.7825; 0.0284; 0.2172; 0.1444; -1.1779];          % get starting th's 
xg = [-0.7825; 0.0284; 0.2172; 0.1444; -1.1779];       % get end th's

xR=[];
xR(:,1)=[x0; zeros(5,1)];
horizon = 30;                       % specify horizon/ time duration
% Generatie reference
% line in the state space
x_ = [];
for i=1:nstate/2
    x_ = [x_; linspace(x0(i),xg(i),horizon+1)];
end%
W = zeros(5,horizon+1);             % penalize angular velocity
x_ = [x_;W];
xref_ = [];
for i = 1:horizon+1%
    xref_ = [xref_; x_(:,i)];
end   
xori=xref_(nstate+1:end);
x_ = xori;
% zeros input
uori = zeros(horizon*njoint,1);

for cc = 1
%% Environment
% Set up obstacles
obs{1}.shape = 'cylinder';
obs{1}.A = [1 0 ; 0 1];
obs{1}.l = [[3806;8413;1]./1000 [3606;8413;1038]./1000]; % cylinder (straight line);% % for M200i
obs{1}.D = 0.2;
obs{1}.epsilon =  0.25;

%% Cost Fn Parameters
tic
Aaug=[];Baug=zeros(horizon*nstate,horizon*nu);Qaug=zeros(horizon*nstate);
Q=[];
Q(1:njoint,1:njoint)=[10 0 0 0 0;
    0 10 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];
Q(1:njoint,njoint+1:2*njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,1:njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,njoint+1:2*njoint)=[10 0 0 0 0;
    0 10 0 0 0;
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
QQ = Baug'*Qaug*Baug+R.*50;
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
sys_info.MAX_input = kron(ones(horizon,1),[1;1;pi;pi;pi]*robot.delta_t);

%% Solver
% baseline
eval = EVAL(sys_info);
Cost_b = eval.get_Cost_b();

for i=1:size(SOLVER,2)
    % select solver
    if ~isempty(SOLVER{i})
        switch(SOLVER{i}.name)            
            case 'PSGCFS'
                % PSGCFS
                tic
                self = PSGCFS_FANUC(obs,sys_info,ROBOT);
                self = self.optimizer();
                Time_{2,cc} = toc; 
                Cost{2,cc} = self.eval.cost_all;
                e_u_all{2,cc} = self.eval.e_u_all;
                e_cost_all{2,cc} = self.eval.e_cost_all;
            case 'CFS'
                % CFS
                tic
                self = CFS_FANUC(obs,sys_info,ROBOT);
                self = self.optimizer();
                Time_{3,cc} = toc; 
                Cost{3,cc} = self.eval.cost_all;
                e_u_all{3,cc} = self.eval.e_u_all;
                e_cost_all{3,cc} = self.eval.e_cost_all;
                cost_all = self.eval.cost_all;
        end

        % get data
        uref = self.u;
        x_ = self.x_;
        xref_ = reshape(x_,[nstate,horizon]);        
    end
end
Time_{1} =toc;

%% Visualization 
% robot animation
xref_ = reshape(x_,[10,horizon]);
figure(1)
[all_ee] = plot_manipulator2(ROBOT,obs,xref_,robot,njoint);

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
        '\theta_{1,ref}','\theta_{2,ref}','\theta_{3,ref}','\theta_{4,ref}','\theta_{5,ref}','Location','northoutside','Orientation','horizontal','NumColumns',5);
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Input [rad]')

%% plot cost all
figure(4)
plot(self.eval.cost_all)
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Cost')

end