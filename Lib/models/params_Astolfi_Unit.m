%% params_pendulum
% create params structure for pendulum 
function params = params_Astolfi_Unit

    % pendulum parameters
    params.A = [0 1; 0 0];
    params.B = [0; 1];
    params.C = [1, 0];
   
    % number of reference trajectories (>1 for control design)
    params.Ntraj = 1;
    
    % reference init
    params.X(1).val(:,1) = [2; 3];

    for traj=2:params.Ntraj
        params.X(traj).val(:,1) = params.X(traj-1).val(:,1);
    end
    
    % position in the state vector of the parameters
    params.estimated_params = [];
    
    % which vars am I optimising
    params.opt_vars = [1:2];
    
    % not opt vars
    tmp = 1:length(params.X(1).val(:,1));
    tmp_idx = tmp;
    for i=1:length(params.opt_vars)
        tmp_idx = intersect(tmp_idx,find(tmp~=params.opt_vars(i)));
    end
    params.nonopt_vars = tmp_idx;
end