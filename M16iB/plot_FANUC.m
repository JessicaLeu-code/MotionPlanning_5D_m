function [all_ee] = plot_FANUC(obs,xref,robot,njoint,D,horizon,offset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate motion planning animetion for FANUC M16iB
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs = obs.*1000;
r = D*700;
[X,Y,Z] = cylinder(r);
surf(X+obs(1,1),Y+obs(2,1),Z*(obs(3,2)-obs(3,1))+obs(3,1));
hold on
[X,Y,Z] = sphere;
surf(X*r+obs(1,1),Y*r+obs(2,1),Z*r+obs(3,1));
surf(X*r+obs(1,1),Y*r+obs(2,1),Z*r+obs(3,2));

plot3(obs(1,1),obs(2,1),obs(3,1),'*')
plot3(obs(1,2),obs(2,2),obs(3,2),'*')
% pause
h=[];

all_ee = [];
for i=1:1:horizon
    pause(0.05)
    delete(h);
    robot.DH(1:njoint,1)=xref(1:5,i);
        
    % Use the following lines to draw capsules
%     color=[i/horizon,i/horizon,i/horizon];
%     valpha=0.2;
%     RobotCapLink;
    
    % Use the following lines to draw robot arm
    valpha=1;%i/horizon;
    RobotFigureLink;
    all_ee = [all_ee  v(end,:)'];
end

end