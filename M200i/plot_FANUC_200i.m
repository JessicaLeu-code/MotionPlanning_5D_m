function [all_ee] = plot_FANUC_200i(obs,xref,robot,njoint,horizon,offset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate motion planning animetion for FANUC M200i
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    RobotFigureLink_200i;
    all_ee = [all_ee  v(end,:)'];
end

end