ROBOT = 'LRMate200iD7L';
load(['figure/',ROBOT,'.mat']);
offset= [3150, 8500,330]./1000;
factor = 1;
%% links
M={};
M{1}=[eye(3) offset';zeros(1,3) 1];

DH=robot.DH;
%DH(2,1) = -(DH(2,1)-pi/2);
DH(2,1) = DH(2,1)- pi/2;
%DH(4,1) = -DH(4,1);
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


