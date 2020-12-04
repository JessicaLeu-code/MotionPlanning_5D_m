function [all_ee] = plot_200i(obs,xref,robot,njoint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate motion planning animetion for FANUC M200i
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:size(obs,2)
    if ~isempty(obs{j})
        switch obs{j}.shape
            case 'cylinder'            
                r = obs{j}.D/2;
                [X,Y,Z] = cylinder(r);
                surf(X+obs{j}.l(1,1),Y+obs{j}.l(2,1),Z*(obs{j}.l(3,2)-obs{j}.l(3,1))+obs{j}.l(3,1));
                hold on
                [X,Y,Z] = sphere;
                surf(X*r+obs{j}.l(1,1),Y*r+obs{j}.l(2,1),Z*r+obs{j}.l(3,1));
                surf(X*r+obs{j}.l(1,1),Y*r+obs{j}.l(2,1),Z*r+obs{j}.l(3,2));
                
            case 'rectangle'            
                ob = Polyhedron('V',obs{j}.poly');
                plot_marker(j) = ob.plot('color','g');
                hold on
        end
    end    
end

% pause
h=[];

all_ee = [];
horizon = size(xref,2);
for j=1:1:horizon
    pause(0.05)
    delete(h);
    robot.DH(1:njoint,1)=xref(1:5,j);
        
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