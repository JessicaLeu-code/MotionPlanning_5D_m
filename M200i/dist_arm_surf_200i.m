function [d, linkid] = dist_arm_3D_200i_2(theta,base,obs,robot)


nstate = size(theta,1);
DH = robot.DH(1:nstate,:);
for i=1:nstate
DH(i,1)=theta(i);
end%theta,d,a,alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%% for 200i
DH(2,1) = DH(2,1)- pi/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = Inf;
if size(base,2)>1
    base=base';
end
cap = robot.cap;
pos=CapPos(base,DH,cap);

for i=1:nstate
    [dis, points] = point2surface_dis(pos{i}.p,obs);
    if norm(dis)<0.0001
        dis = -norm(points(:,1)-pos{i}.p(:,2));
    end        
    if dis < d
        d = dis;
        linkid=i;
    end
end
end