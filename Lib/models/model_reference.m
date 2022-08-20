%% MODEL_REFERENCE
% file: model_reference.m
% author: Federico Oliva
% date: 22/06/2022
% description: this function describes the dynamics equation to be used as
% reference in the control design
% INPUT:
% t: time instant
% x: state vector
% params: structure with all the necessary parameters 
% obs: observer class instance (may be not used)
% OUTPUT:
% x_dot: dynamics equations
function x_dot = model_reference(t, x, params, obs)

    % init the dynamics
    x_dot = zeros(length(x),1);    
    
    % compute the control
    tdiff = obs.setup.time-t;   
    pos = find(abs(tdiff) == min(abs(tdiff)),1,'first');    
    % check length of measurements
    pos = max(1,min(pos,size(obs.init.Y_full_story(obs.init.traj).val,3)));
    drive_out = drive(obs,obs.init.X(obs.init.traj).val(:,pos),x,obs.init.Y_full_story(obs.init.traj).val(:,:,max(1,pos-1):pos),pos);
    params.u = params.input(t,drive_out,params,obs);   
    
    % save input
    obs.init.input_story(obs.init.traj).val(:,pos) = params.u;
    
    % Plant - true one
    Ap = [params.A1 params.A2; params.A3 params.A4];
    Bp = [params.B1; params.B2]; 
         
    % Control
    Ac = [0 1; params.a0 params.a1];
    Bc = [params.b0; params.b1];
    
    % u1 = uc (input to plant - ref track)
    % u2 = ec (input to controller - ref track)
    % u3 = ur (input to reference model - ref track)
    % u4 = ui (input to plant - sys id)
    
    %%% model dynamics %%%           
    % plant dynamics - reference tracking
    x_dot(1:2,:) = Ap*x(1:2,:) + Bp*(params.u(1,:));
    % controller dynamics - reference tracking
    x_dot(3:4,:) = Ac*x(3:4,:) + Bc*(params.u(2,:));
    % reference model dynamics - reference tracking (useless here)
    x_dot(5,:) = params.alpha*x(5,:)+abs(params.alpha)*params.u(3,:);
    % plant dynamics - system identification
    x_dot(6:7,:) = Ap*x(6:7,:) + Bp*(params.u(4,:));
        
end