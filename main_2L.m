%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planning comparison on manipulator
% With planning algorithm: CFS, PSGCFS
% Robot model: Two-link robot
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

%% Load robot info
ROBOT = '2L';                   % 'M200i', '2L'
robot = robotproperty2(ROBOT);

njoint=2;                          % first five joint angles
nstate=4;                         % [ th x 5 ; w x 5 ]
nu=2;                              % [ alpha x 5 ]
DH=robot.DH;

%% Solver set up
SOLVER{1}.name ='CFS';
SOLVER{1}.plot = true;
% SOLVER{2}.name = 'PSGCFS';
% SOLVER{2}.plot = true;


%% Initialization  
% 2L
x0 = [0; 0];              % get starting th's 
xg = [pi/2; 0];       % get end th's

xR=[];
xR(:,1)=[x0; zeros(njoint,1)];
horizon = 40;                       % specify horizon/ time duration
% Generatie reference
% line in the state space
x_ = [];
for i=1:nstate/2
    x_ = [x_; linspace(x0(i),xg(i),horizon+1)];
end%
W = zeros(njoint,horizon+1);             % penalize angular velocity
x_ = [x_;W];
xref_ = [];
for i = 1:horizon+1%
    xref_ = [xref_; x_(:,i)];
end   
xori=xref_(nstate+1:end);
%x_ = xori;
x_ = kron(ones(horizon,1),[x0;0;0]);     % stationary reference
% zeros input
uori = zeros(horizon*njoint,1);

Time_ = {};
for cc = 1
%% Environment
% Set up obstacles
obs{1}.shape = 'circle';
obs{1}.A = [1 0 ; 0 1];
obs{1}.c = [0.3;0.3];
ccc = obs{1}.c;
obs{1}.l = [[ccc(1);ccc(2);0] [ccc(1);ccc(2);0]];% % for M200i
obs{1}.D = 0.05;
obs{1}.epsilon =  0.05;

%obs{2}.l = [[4106;8313;1]./1000 [4106;8313;1438]./1000];%
view_area = [-0.5  0.5 -0.5 0.5];
%% Cost Fn Parameters
Aaug=[];Baug=zeros(horizon*nstate,horizon*nu);Qaug=zeros(horizon*nstate);
Q=[];
Q(1:njoint,1:njoint)=[10 0;
                       0 1 ];
Q(1:njoint,njoint+1:2*njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,1:njoint)=0.1*eye(njoint);
Q(njoint+1:2*njoint,njoint+1:2*njoint)=[10 0;
                                        0 1 ];
for i=1:horizon
    Aaug=[Aaug;robot.A([1:njoint,njoint+2:1+2*njoint],[1:njoint,njoint+2:1+2*njoint])^i];
    Qaug((i-1)*nstate+1:i*nstate,(i-1)*nstate+1:i*nstate)=Q*0.1;
    if i==horizon
        Qaug((i-1)*nstate+1:i*nstate,(i-1)*nstate+1:i*nstate)=Q*10000;
    end
    for j=1:i
        Baug((i-1)*nstate+1:i*nstate,(j-1)*nu+1:j*nu)=robot.A([1:njoint,njoint+2:1+2*njoint],[1:njoint,njoint+2:1+2*njoint])^(i-j)*robot.B([1:njoint,njoint+2:1+2*njoint],1:nu);
    end
end
R=eye(horizon*nu);
for i=1:horizon
    R((i-1)*nu+1:i*nu,(i-1)*nu+1:i*nu)=[5 0 ;
                                        0 4];
end
R=R+R';
QQ = Baug'*Qaug*Baug+R.*0.1;
gaug = kron(ones(horizon,1),[xg;0;0]);
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
sys_info.x_ = x_;
sys_info.alpha = 1/max(svd(QQ));
sys_info.lim = [0.1;0.2];
% Solver setup
sys_info.options = optimoptions('quadprog','Display','off');
sys_info.alpha = 1/max(svd(QQ)); 
sys_info.epsilon_O = 1e-6;
sys_info.MAX_O_ITER = 100;  % 6 for chomp
sys_info.MAX_input = kron(ones(horizon,1),[1;1]*0.5*robot.delta_t);

%% Solver
% baseline
eval = EVAL(sys_info);
Cost_b = eval.get_Cost_b();

%%%%%%%%%%%%%%%
s_Solver;
%%%%%%%%%%%%%%%

%% Visualization 
% plot inputs
uref_ = reshape(uref,[njoint,horizon]);
uori_ = reshape(uori,[njoint,horizon]);
figure(2);
plot(uref_','--')
legend('\alpha_1','\alpha_2')
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Input [rad/s^2]')

% plot positions
x_ = reshape(x_,[nstate,horizon]);
figure(3)
plot(x_(1:njoint,:)')
hold on
legend('1','2','3','4','5')
% plot ref
xori_ = reshape(xori,[nstate,horizon]);
plot(xori_(1:njoint,:)','--')
legend('\theta_1','\theta_2','\theta_{1,ref}','\theta_{2,ref}',...
    'Location','northoutside','Orientation','horizontal','NumColumns',5);
xlabel(['Time steps ' '[' num2str(robot.delta_t) 's/step]'])
ylabel('Input [rad]')


end