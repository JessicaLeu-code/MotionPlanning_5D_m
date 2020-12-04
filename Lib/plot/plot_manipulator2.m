function [all_ee] = plot_manipulator2(ROBOT,obs,route,robot,njoint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot experiment result (robot and obstacles)
% Robot model: FANUC M200i, M16iB, 2-Link robot
% ROBOT: robot name,                            (char)'name'
% obs: obstacle info,                           (cell) obs{}.n
% route: joint angle planning result,           (vector)[#joints x horizon]
% robot: robot info,                            (struct)'robot.n'
% njoint: #joints involves in planning,         (int) n
%
% Jessica Leu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
horizon = size(route,2);
switch ROBOT
    case 'M200i'        
        DrawMap2;
        delete(h)
        [all_ee] = plot_200i(obs,route,robot,njoint);       
    
    case 'M16iB'
        DrawMap;
        delete(h)
        [all_ee] = plot_M16iB(obs,route,robot,njoint);
        
    case '2L'
        view_area = [-0.5  0.5 -0.5 0.5];
        [all_ee] = plot_2L(obs,route,robot,njoint,view_area);        
end
end
