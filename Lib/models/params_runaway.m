%% runaway init model
function params = params_runaway

    % plant data

    % electric charge
    params.Q = 15;

    % spontaneous emission
    params.S = 0;

    % parameters
    params.gamma = 2;
    params.gamma1 = 10;
    params.ni = 20;
    params.Wt = 0.3;

    % eps_coef
    params.eps_coef = 1;

    % initial condition
    params.T0 = 5;
    params.W0 = 1;
    
    % initial condition
    params.X = [params.T0; params.W0];
    
    % position in the state vector of the parameters
    params.estimated_params = [];
end