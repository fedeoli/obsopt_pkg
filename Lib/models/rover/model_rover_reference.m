%% MODEL_ROVER
% file: model_rover.m
% author: Federico Oliva
% date: 10/01/2022
% description: this function describes the dynamics equation of a rover
% t: time instant
% x: state vector
% params: structure with all the necessary parameters 
% obs: observer class instance (may be not used)
% OUTPUT:
% x_dot: dynamics equations
function x_dot = model_rover_reference(tspan,x,params,obs)

    % init the dynamics 
    x_dot = zeros(length(x),1);

    % compute the time index
    for i=1:length(tspan)
        tdiff = obs.setup.time-tspan(i);   
        pos(i) = find(abs(tdiff) == min(abs(tdiff)),1,'first');    
        pos(i) = max(1,pos(i));        
    end    
    
    % compute the control
    params.u = params.input(tspan,x,params);    
    obs.init.input_story_ref(obs.init.traj).val(:,pos(1)) = params.u(:,1);    
    
    % model dynamics
    % x axis
    x_dot(1) = x(2);
    x_dot(2) = x(3) + params.u(1,1);    
    
    % y axis
    x_dot(6) = x(7);
    x_dot(7) = x(8) + params.u(2,1);    

    % all the remaining are the anchors
    
end