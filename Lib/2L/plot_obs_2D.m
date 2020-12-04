function [plot_marker] = plot_obs_2D(obs,view_area)
num_obs = size(obs,2); 
for i=1:num_obs
    if ~isempty(obs{i})
        switch obs{i}.shape
            case 'circle'            
                [x, y, z] = ellipsoid(obs{i}.c(1), obs{i}.c(2),0,obs{i}.m/sqrt(obs{i}.A(1,1)),obs{i}.m/sqrt(obs{i}.A(2,2)),0,30);
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

end
            
            
            
            