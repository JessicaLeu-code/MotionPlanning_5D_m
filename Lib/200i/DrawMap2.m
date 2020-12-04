%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Plot for experiment result animation (robot and obstacles)
% Robot model: FANUC M200i
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
offset= [3150, 8500,330]./1000;
%load('figure/stage.mat');
%patch('Faces',stage.f,'Vertices',stage.v,'FaceVertexCData',stage.c,'FaceColor',stage.color,'EdgeColor','None');

ROBOT1 = 'LRMate200iD7L';
load(['figure/',ROBOT1,'.mat']);
% % robot=robotproperty2(ROBOT);
% % hold on
% % offset= robot.base';



factor = 1;
% factor = 1/1000; % For M16iB 
h=[];
%% Base
for i=1:1
    f=base{i}.f; v=base{i}.v*factor+ones(size(base{i}.v,1),1)*offset; color=base{i}.color;
    patch('Faces',f,'Vertices',v,'FaceColor',color,'EdgeColor','None');
end

%% links
M={};
M{1}=[eye(3) offset';zeros(1,3) 1];
DH=robot.DH;
%DH(:,1) = xref(1:njoint);
DH(:,1) = [0;0;0;0;0;pi/3];
DH(2,1) = -(DH(2,1)-pi/2);
DH(2,1) = DH(2,1)- pi/2;
DH(4,1) = -DH(4,1);
for i=1:6
    R=[cos(DH(i,1)) -sin(DH(i,1))*cos(DH(i,4)) sin(DH(i,1))*sin(DH(i,4));
        sin(DH(i,1)) cos(DH(i,1))*cos(DH(i,4)) -cos(DH(i,1))*sin(DH(i,4));
        0  sin(DH(i,4)) cos(DH(i,4))];
    T=[DH(i,3)*cos(DH(i,1));DH(i,3)*sin(DH(i,1));DH(i,2)];
    M{i+1}=M{i}*[R T;zeros(1,3) 1];
end
%%
for i=1:5
    v=link{i}.v*factor; f=link{i}.f;  color=link{i}.color;
    for j=1:size(v,1)
        v(j,:)=v(j,:)*M{i+1}(1:3,1:3)'+M{i+1}(1:3,4)';
    end
    h(i)=patch('Faces',f,'Vertices',v,'FaceColor',color,'EdgeColor','None');
    
end

%% Payload
i=6;%c=payload{1}.c
v=payload{1}.v;f=payload{1}.f;;color=payload{1}.color;
for j=1:size(v,1)
    v(j,:)=v(j,:)*M{i+1}(1:3,1:3)'+M{i+1}(1:3,4)';
end
h(i)=patch('Faces',f,'Vertices',v,'FaceColor',color,'EdgeColor','None');


%%
% axis equal
% % axis([-1 1,-1 1,0,2])
% view([1,0.4,1])
% %lighting flat
% light=camlight('right');
axis equal
xlim = [-0.4, 1.2]+offset(1);
ylim = [-1.2, 1]+offset(2);
zlim = [-0, 1.5];%+offset(3);

axis([xlim, ylim, zlim])
% view([1,0.4,1])
view([1,-0.5,0.4])
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]'); 
lighting flat
light=camlight('left');


wall{1}.handle=fill3([xlim(1),xlim(1),xlim(2),xlim(2)],[ylim(1),ylim(2),ylim(2),ylim(1)],[zlim(1),zlim(1),zlim(1),zlim(1)],[0.8 0.8 0.8]);
wall{2}.handle=fill3([xlim(1),xlim(1),xlim(1),xlim(1)],[ylim(1),ylim(1),ylim(2),ylim(2)],[zlim(1),zlim(2),zlim(2),zlim(1)],[0,0.9,0.9]);
% wall{3}.handle=fill3([xlim(1),xlim(1),xlim(2),xlim(2)],[ylim(2),ylim(2),ylim(2),ylim(2)],[zlim(1),zlim(2),zlim(2),zlim(1)],[0,0.9,0.9]);
%%
% load('figure/block.mat');
% for i=1:1
%     %c=block{i}.c;
%     f=block{i}.f; v=block{i}.v*factor+ones(size(block{i}.v,1),1)*[0.45;0.24;0.29]';  color=block{i}.color;
%     patch('Faces',f,'Vertices',v,'FaceColor',color,'EdgeColor','None');
% end

% %% Obs
% D = 0.2;
% 
% r = D*700;
% %% Plot enclosed cylinder
% 
% % Sample values
% %obs = obs.*1000;
% h1 = (obs(3,2)-obs(3,1))+obs(3,1);     % height
% ra = D/2;   % radius
% % Create constant vectors
% tht = linspace(0,2*pi,100); z = linspace(0,h1,20);
% % Create cylinder
% xa = repmat(ra*cos(tht),20,1); ya = repmat(ra*sin(tht),20,1);
% za = repmat(z',1,100);
% % To close the ends
% X = [xa*0; flipud(xa); (xa(1,:))*0]; Y = [ya*0; flipud(ya); (ya(1,:))*0];
% Z = [za; flipud(za); za(1,:)];
% % Draw cylinder
% [TRI,v]= surf2patch(X+obs(1,1),Y+obs(2,1),Z,'triangle'); 
% patch('Vertices',v,'Faces',TRI,'facecolor',[0.5 0.8 0.8],'facealpha',0.8);
%  grid on; axis square; title('Cylinder','FontSize',12)
% 
% % pause
