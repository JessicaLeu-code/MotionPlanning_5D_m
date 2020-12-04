function [cap_all] = get_cap_all(theta_all,base,robot)

horizon = size(theta_all,2);
cap = robot.cap;
nstate = size(theta_all,1);
DH = robot.DH(1:nstate,:);
cap_all = cell(horizon);
for h = 1:horizon
    for i=1:nstate
        DH(i,1)=theta_all(i,h);
    end%theta,d,a,alpha

    
    if size(base,2)>1
        base=base';
    end

    cap_all{h}=CapPos(base,DH,cap);
end

end