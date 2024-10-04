%% SIMULATION_REALDATA_BATTERY_SIMULINK
% file: simulation_general_v3.m
% author: Federico Oliva
% date: 27/02/2023
% description: function to setup and use the MHE observer on battery from
% real data of CALCE dataset 
% INPUT: none
% OUTPUT: params,obs
function [params,obs] = simulation_realdata_battery_calce

% generate from simulink
params_sim = params_battery_simulink_calce;

% init observer buffer (see https://doi.org/10.48550/arXiv.2204.09359)
Nw = 30;
Nts = 30;

% Tushar setup
% Nw = 30; % observer window lenth
% Nts =2; % downsampling factor
% multirate downsampling
% temp_Nts = ones(1,Nw) *2;
% temp_Nts(1:Nw-5)=20; 
% Nts = temp_Nts;

% noise
rng default

% set sampling time
Ts = params_sim.Ts;

% set initial and final time instant
% remember to set the final time and sampling time accordingly to the data
% that you measured
t0 = params_sim.time(1);
tend = params_sim.time(end);
% uncomment to test the MHE with a single optimisation step
% tend = 5*(Nw*Nts-1)*Ts;

%%%%%% general functions
params_init = @params_battery_calce;
params_update = @params_update_battery_tushar;
model = @model_battery_calce;
% model_reference = @model_battery_calce;
measure = @measure_battery_tushar;
% measure_reference = @measure_battery_tushar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% filters %%%%
[filter, filterScale, reference] = filter_define(Ts,Nts);

%%%% integration method %%%%
% ode45-like integration method. For discrete time systems use @odeDD
ode = @odeDD;

%%%% input law %%%
input_law = @control_battery;

%%%% params init %%%%
params = model_init('Ts',Ts,'T0',[t0, tend],'noise',0, 'params_update', params_update, ...
            'model',model,'measure',measure,'ode',ode, 'odeset', [1e-3 1e-6], ...
            'input_enable',1,'input_law',input_law,'params_init',params_init,'addons',params_sim);
             
%%%% observer init %%%%
%%%% define arrival cost %%%%%
terminal_states = params.opt_vars;
terminal_weights = 1e0*ones(size(terminal_states));
% SOC%
terminal_weights(1) = 5e1;
% OCV %
terminal_weights([3 6]) = 1;
% R0 %
% terminal_weights([4 8 12]) = 1e-2;
% OCV %
% terminal_weights([5 9 13]) = 1e-2;
% OCV %
% terminal_weights([6 10 14]) = 1e4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% create observer class instance. For more information on the setup
% options check directly the class constructor in obsopt.m
obs = obsopt('DataType', 'real', 'optimise', 1, 'MultiStart', params.multistart, 'J_normalise', 1, 'MaxOptTime', Inf, ... 
          'Nw', Nw, 'Nts', Nts, 'ode', ode, 'PE_maxiter', 0, 'WaitAllBuffer', 0, 'params',params, 'filters', filterScale,'filterTF', filter, ...
          'Jdot_thresh',0.95,'MaxIter', 2, 'Jterm_store', 1, 'AlwaysOpt', 1 , 'print', 0 , 'SafetyDensity', Inf, 'AdaptiveParams', [10 20 2 20 1 0 0 params.OutDim_compare], ...
          'AdaptiveSampling',0, 'FlushBuffer', 1, 'opt', @fminsearchcon, 'terminal', 1, 'terminal_states', terminal_states, 'terminal_weights', terminal_weights, 'terminal_normalise', 1, ...
          'ConPos', [], 'LBcon', [], 'UBcon', [],'NONCOLcon',@nonlcon_fcn,'Bounds', 1,'BoundsPos',[1 4 5],'BoundsValLow',[1e-3 1e-3 1e-3],'BoundsValUp',[1 1e3 1e3]);



%% %%%% SIMULATION %%%%
% generate data
%%%% SWITCH Y WITH THE COLLECTED DATA %%%%%%
Niter = obs.setup.Niter;
Ntraj = obs.setup.Ntraj;

% start time counter
t0 = tic;

% integration loop
for i = 1:Niter
    
    % Display iteration step
    if ((mod(i,100) == 0) || (i == 1))
        clc
        disp('MHE:')
        disp(['Iteration Number: ', num2str(obs.setup.time(i)),'/',num2str(obs.setup.time(Niter))])
        disp(['Last J: ',num2str(obs.init.Jstory(end))]);        
    end
    
    % set current iteration in the obsopt class
    obs.init.ActualTimeIndex = i;
    obs.init.t = obs.setup.time(obs.init.ActualTimeIndex);    
    
    %%%% PROPAGATION %%%%
    % forward propagation of the previous estimate
    for traj = 1:Ntraj
        
        % set trajectories
        obs.init.traj = traj;
        obs_slow.init.traj = traj;
                 
        % propagate only if the time gets over the initial time instant
        if (obs.init.ActualTimeIndex > 1)
            
            % define the time span of the integration
            startpos = obs.init.ActualTimeIndex-1;
            stoppos = obs.init.ActualTimeIndex;
            tspan = obs.setup.time(startpos:stoppos); 
            
            % fast
            X = obs.setup.ode(@(t,x)obs.setup.model(t, x, obs.init.params, obs), tspan, obs.init.X_est(traj).val(:,startpos),params.odeset);
            obs.init.X_est(traj).val(:,startpos:stoppos) = [X.y(:,1),X.y(:,end)];
                                                      
        end                
        
        %%%% FILTER MEASUREMENT %%%%                 
        % filter on the measurements    
        % reference filter (filterScale)
        y_meas(traj).val = params.y_sim(:,i);            
    end
       
    %%%% OBSERVER %%%%
    
    % call the observer                 
    tfast = tic;
    obs = obs.observer(obs.init.X_est,y_meas);            
    obs.init.iter_time(obs.init.ActualTimeIndex) = toc(tfast);                          

end
% 
% overall computation time
obs.init.total_time = toc(t0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% post sim operations %%%
obs.init.X.val = zeros(params.dim_state,params.Niter);
obs.init.X.val(1,:) = params_sim.out.simout.ECM_soc.Data';
obs.init.X.val(3,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.OCV,params_sim.out.simout.ECM_soc.Data);
obs.init.X.val(4,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.R0,params_sim.out.simout.ECM_soc.Data);
obs.init.X.val(5,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.R1,params_sim.out.simout.ECM_soc.Data);
obs.init.X.val(6,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.C1,params_sim.out.simout.ECM_soc.Data);

obs.init.Ytrue_full_story.val = zeros(obs.setup.Nfilt,params.OutDim,params.Niter);
obs.init.Ytrue_full_story.val(1,1,:) = params_sim.out.simout.ECM_Vb.Data';
% the whole process could be long, why not going for a nap? No worries, 
% this "sounds" like a nice way to wake up. (Uncomment)
% load handel
% sound(y,Fs)
obs.init.Nw_Nts=Nts*Nw;
end

