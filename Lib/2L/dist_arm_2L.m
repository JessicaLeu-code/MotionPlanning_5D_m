function [d, linkid] = dist_arm_2L(theta,base,obs,robot)


nstate = size(theta,1);

d = Inf;
if size(base,2)>1
    base=base';
end

[ pos] = CapPos2(theta,base,robot);

for i=1:nstate
    [dis, points] = distLinSeg(pos{i}.p(:,1),pos{i}.p(:,2), obs(:,1),obs(:,2));
    if norm(dis)<0.0001
        dis = -norm(points(1:3,1)-pos{i}.p(:,2));
    end        
    if dis < d
        d = dis;
        linkid=i;
    end
end
end