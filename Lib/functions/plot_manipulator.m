function [all_ee] = plot_manipulator(ROBOT,obs,route,robot,njoint,D,path_length,offset)
switch ROBOT
    case 'M200i'
        DrawMap2;
        delete(h)
        [all_ee] = plot_200i(obs,route,robot,njoint);
        
    
    case 'M16iB'
        horizon = path_length;
        DrawMap;
        delete(h)
        [all_ee] = plot_FANUC(obs{1}.l,route,robot,njoint,D,path_length,offset);
    case '2L'
        horizon = path_length;
        DrawMap;
        delete(h)
        [all_ee] = plot_2L(obs,route,robot,njoint,view_area);        
end
end
