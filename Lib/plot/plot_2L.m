function [end_all] = plot_2L(obs,xref,robot,njoint,view_area)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate motion planning animetion for Two-link robot
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot obs
plot_marker=[];
num_obs = size(obs,2); 
for i=1:num_obs
    if ~isempty(obs{i})
        switch obs{i}.shape
            case 'circle'            
                [x, y, z] = ellipsoid(obs{i}.c(1), obs{i}.c(2),0,obs{i}.D/sqrt(obs{i}.A(1,1)),obs{i}.D/sqrt(obs{i}.A(2,2)),0,30);
                plot_marker(i) = surf(x, y, z);
                hold on
            case 'rectangle'            
                ob = Polyhedron('V',obs{i}.poly');
                plot_marker(i) = ob.plot('color','g');
                hold on
        end
    end
    axis(view_area)
    axis equal
    view(2)
end

DH=robot.DH(1:njoint,:);

end_all = [];
h= [];
links_ =[];
etrj = [];
cons_c = 1;
horizon = size(xref,2);
for i=1:1:horizon
    pause(0.1)
    delete(h);
    delete(links_);
    delete(etrj);
    theta=xref(1:njoint,i);
    for j=1:njoint
        DH(j,1)=theta(j);
    end%theta,d,a,alpha
    [pos,M]=CapPos2(theta,robot.base,robot);
    
    %% plot arm
  
    for j = 1:njoint
        if j<njoint
       links_(j) = plot3([ pos{j}.p(1,1), pos{j+1}.p(1,1)],[ pos{j}.p(2,1), pos{j+1}.p(2,1)],[ pos{j}.p(3,1), pos{j+1}.p(3,1)],'k-','color',[1-cons_c,1-cons_c,1-cons_c/1.5],'LineWidth',3);
        end  
       
       etrj(j) = plot3( pos{j}.p(1,2), pos{j}.p(2,2), pos{j}.p(3,2),'o-','color',[1-cons_c/3.5,1-cons_c/2.5,1-cons_c],'LineWidth',3);
       hold on
    end
    links_(j) = plot3([ pos{njoint}.p(1,1), pos{njoint}.p(1,2)],[ pos{njoint}.p(2,1), pos{njoint}.p(2,2)],[ pos{njoint}.p(3,1), pos{njoint}.p(3,2)],'k-','color',[1-cons_c,1-cons_c,1-cons_c/1.5],'LineWidth',3);
    end_all = [end_all pos{j}.p(1:3,2)];    
    xlabel(['x[m] Time step:' num2str(i)])
    ylabel('y[m]')
    zlabel('z[m]')
    axis equal
    axis([view_area])
    view(2)
end

end