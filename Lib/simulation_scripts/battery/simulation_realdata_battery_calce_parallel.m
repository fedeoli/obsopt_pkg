%% SIMULATION_REALDATA_BATTERY_SIMULINK
% file: simulation_general_v3.m
% author: Federico Oliva
% date: 27/02/2023
% description: function to setup and use the MHE observer on battery from
% real data of CALCE dataset 
% INPUT: none
% OUTPUT: params,obs
function [params_out,obs] = simulation_realdata_battery_calce_parallel

% generate from simulink
params_fast = params_battery_simulink_calce;
params_sim = params_fast;
clear params_fast

% init observer buffer (see https://doi.org/10.48550/arXiv.2204.09359)
% Nw = 30;

% noise
rng default

% set sampling time
Ts = params_sim.Ts;

% set initial and final time instant
% remember to set the final time and sampling time accordingly to the data
% that you measured
t0 = params_sim.time(1);
tend = params_sim.time(end);
% tend = 1000;

%%%%%% general functions
params_update = @params_update_battery_tushar;
model = @model_battery_calce;
% model_reference = @model_battery_calce;
measure = @measure_battery_tushar;
% measure_reference = @measure_battery_tushar;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%% integration method %%%%
% ode45-like integration method. For discrete time systems use @odeDD
ode = @odeDD;

%%%% input law %%%
input_law = @control_battery;

%%% FAST OBSOPT %%%

% params init - fast
params_fast = model_init('Ts',Ts,'T0',[t0, tend],'noise',0, 'params_update', params_update, ...
            'model',model,'measure',measure,'ode',ode, 'odeset', [1e-3 1e-6], ...
            'input_enable',1,'input_law',input_law,'params_init',@params_battery_calce_fast,'addons',params_sim);             

% create observer class instance. For more information on the setup
% options check directly the class constructor in obsopt.m

% instance for SOC, R0, R1
Nw_fast = 30;
% NTs
Nts_fast = 1;

% filters
[filter_fast, filterScale_fast, ~] = filter_define(Ts,Nts_fast);
% terminal states
terminal_states_fast = params_fast.opt_vars;
terminal_weights_fast = 1e0*ones(size(terminal_states_fast));

% class instance
obs_fast = obsopt('DataType', 'real', 'optimise', 1, 'MultiStart', params_fast.multistart, 'J_normalise', 1, 'MaxOptTime', Inf, ... 
          'Nw', Nw_fast, 'Nts', Nts_fast, 'ode', ode, 'PE_maxiter', 0, 'WaitAllBuffer', 0, 'params',params_fast, 'filters', filterScale_fast,'filterTF', filter_fast, ...
          'Jdot_thresh',0.95,'MaxIter', 2, 'Jterm_store', 1, 'AlwaysOpt', 1 , 'print', 0 , 'SafetyDensity', Inf, 'AdaptiveParams', [10 20 2 5 1 0 0 params_fast.OutDim_compare], ...
          'AdaptiveSampling',1, 'FlushBuffer', 1, 'opt', @fminsearchcon, 'terminal', 1, 'terminal_states', terminal_states_fast, 'terminal_weights', terminal_weights_fast, 'terminal_normalise', 1, ...
          'ConPos', [], 'LBcon', [], 'UBcon', [],'NONCOLcon',@nonlcon_fcn,'Bounds', 1,'BoundsPos',[1 4 5],'BoundsValLow',[1e-3 1e-3 1e-3],'BoundsValUp',[1 1e3 1e3]);

%%% SLOW OBSOPT %%%

% params init - slow
params_slow = model_init('Ts',Ts,'T0',[t0, tend],'noise',0, 'params_update', params_update, ...
            'model',model,'measure',measure,'ode',ode, 'odeset', [1e-3 1e-6], ...
            'input_enable',1,'input_law',input_law,'params_init',@params_battery_calce_slow,'addons',params_sim);             

% create observer class instance. For more information on the setup
% options check directly the class constructor in obsopt.m

% instance for SOC, R0, R1
Nw_slow = 30;
% NTs
Nts_slow = 100;

% Tushar setup
% temp_Nts = ones(1,Nw) *2;
% temp_Nts(1:Nw-5)=20;
% Nts_slow = temp_Nts;

% filters
[filter_slow, filterScale_slow, ~] = filter_define(Ts,Nts_slow);
% terminal states
terminal_states_slow = params_slow.opt_vars;
terminal_weights_slow = 1e0*ones(size(terminal_states_slow));
terminal_weights_slow(1) = 5e1;

% class instance
obs_slow = obsopt('DataType', 'real', 'optimise', 1, 'MultiStart', params_slow.multistart, 'J_normalise', 1, 'MaxOptTime', Inf, ... 
          'Nw', Nw_slow, 'Nts', Nts_slow, 'ode', ode, 'PE_maxiter', 0, 'WaitAllBuffer', 0, 'params',params_slow, 'filters', filterScale_slow,'filterTF', filter_slow, ...
          'Jdot_thresh',0.95,'MaxIter', 2, 'Jterm_store', 1, 'AlwaysOpt', 1 , 'print', 0 , 'SafetyDensity', Inf, 'AdaptiveParams', [10 20 2 50 0.001 0 0 params_slow.OutDim_compare], ...
          'AdaptiveSampling',1, 'FlushBuffer', 1, 'opt', @fminsearchcon, 'terminal', 1, 'terminal_states', terminal_states_slow, 'terminal_weights', terminal_weights_slow, 'terminal_normalise', 1, ...
          'ConPos', [], 'LBcon', [], 'UBcon', [],'NONCOLcon',@nonlcon_fcn,'Bounds', 1,'BoundsPos',[1 4 5],'BoundsValLow',[1e-3 1e-3 1e-3],'BoundsValUp',[1 1e3 1e3]);

%% %%%% SIMULATION %%%%
% generate data
%%%% SWITCH Y WITH THE COLLECTED DATA %%%%%%
Niter = obs_fast.setup.Niter;
Ntraj = obs_fast.setup.Ntraj;

% start time counter
t0 = tic;

% integration loop
for i = 1:Niter
    
    % Display iteration step
    if ((mod(i,100) == 0) || (i == 1))
        clc
        disp('MHE:')
        disp(['Iteration Number: ', num2str(obs_fast.setup.time(i)),'/',num2str(obs_fast.setup.time(Niter))])
        disp(['Last J (fast): ',num2str(obs_fast.init.Jstory(end))]);        
        disp(['Last J (slow): ',num2str(obs_slow.init.Jstory(end))]);  
    end
    
    % set current iteration in the obsopt class
    obs_fast.init.ActualTimeIndex = i;
    obs_fast.init.t = obs_fast.setup.time(obs_fast.init.ActualTimeIndex);    

    obs_slow.init.ActualTimeIndex = i;
    obs_slow.init.t = obs_slow.setup.time(obs_slow.init.ActualTimeIndex);    
    
    %%%% PROPAGATION %%%%
    % forward propagation of the previous estimate
    for traj = 1:Ntraj
        
        % set trajectories
        obs_fast.init.traj = traj;
        obs_fast.init.traj = traj;
                 
        % propagate only if the time gets over the initial time instant
        if (obs_fast.init.ActualTimeIndex > 1)
            
            % define the time span of the integration
            startpos = obs_fast.init.ActualTimeIndex-1;
            stoppos = obs_fast.init.ActualTimeIndex;
            tspan = obs_fast.setup.time(startpos:stoppos); 
            
            % fast
            X = obs_fast.setup.ode(@(t,x)obs_fast.setup.model(t, x, obs_fast.init.params, obs_fast), tspan, obs_fast.init.X_est(traj).val(:,startpos),params_fast.odeset);
            obs_fast.init.X_est(traj).val(:,startpos:stoppos) = [X.y(:,1),X.y(:,end)];

            % slow
            X = obs_slow.setup.ode(@(t,x)obs_slow.setup.model(t, x, obs_slow.init.params, obs_slow), tspan, obs_slow.init.X_est(traj).val(:,startpos),params_slow.odeset);
            obs_slow.init.X_est(traj).val(:,startpos:stoppos) = [X.y(:,1),X.y(:,end)];
                                                      
        end                
        
        %%%% MEASUREMENT %%%%                 
        y_meas(traj).val = params_fast.y_sim(:,i);            
    end
       
    %%%% OBSERVER %%%%
    
    % call the observer                     
    % fast
    tfast = tic;
    obs_fast = obs_fast.observer(obs_fast.init.X_est,y_meas);            
    obs_fast.init.iter_time(obs_fast.init.ActualTimeIndex) = toc(tfast);                          
    % slow
    tslow = tic;
    obs_slow = obs_slow.observer(obs_slow.init.X_est,y_meas);            
    obs_slow.init.iter_time(obs_slow.init.ActualTimeIndex) = toc(tslow);

    % slow updates fast
    if mod(i,Nts_slow) == 0
        obs_fast.init.X_est(traj).val(params_slow.opt_vars,max(1,i-Nw_slow*Nts_slow):i) = obs_slow.init.X_est(traj).val(params_slow.opt_vars,max(1,i-Nw_slow*Nts_slow):i);
        obs_slow.init.X_est(traj).val(params_fast.opt_vars,max(1,i-Nw_fast*Nts_fast):i) = obs_fast.init.X_est(traj).val(params_fast.opt_vars,max(1,i-Nw_fast*Nts_fast):i);
    end

end
% 
% overall computation time
obs_fast.init.total_time = toc(t0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% post sim operations %%%
obs_fast.init.X.val = zeros(params_fast.dim_state,params_fast.Niter);
obs_fast.init.X.val(1,:) = params_sim.out.simout.ECM_soc.Data(1:tend+1)';
obs_fast.init.X.val(3,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.OCV,params_sim.out.simout.ECM_soc.Data(1:tend+1));
obs_fast.init.X.val(4,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.R0,params_sim.out.simout.ECM_soc.Data(1:tend+1));
obs_fast.init.X.val(5,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.R1,params_sim.out.simout.ECM_soc.Data(1:tend+1));
obs_fast.init.X.val(6,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.C1,params_sim.out.simout.ECM_soc.Data(1:tend+1));

obs_fast.init.Ytrue_full_story.val = zeros(obs_fast.setup.Nfilt,params_fast.OutDim,params_fast.Niter);
obs_fast.init.Ytrue_full_story.val(1,1,:) = params_sim.out.simout.ECM_Vb.Data(1:tend+1)';


obs_slow.init.X.val = zeros(params_fast.dim_state,params_fast.Niter);
obs_slow.init.X.val(1,:) = params_sim.out.simout.ECM_soc.Data(1:tend+1)';
obs_slow.init.X.val(3,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.OCV,params_sim.out.simout.ECM_soc.Data(1:tend+1));
obs_slow.init.X.val(4,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.R0,params_sim.out.simout.ECM_soc.Data(1:tend+1));
obs_slow.init.X.val(5,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.R1,params_sim.out.simout.ECM_soc.Data(1:tend+1));
obs_slow.init.X.val(6,:) = interp1(params_sim.input_data.SOC,params_sim.input_data.C1,params_sim.out.simout.ECM_soc.Data(1:tend+1));

obs_slow.init.Ytrue_full_story.val = zeros(obs_fast.setup.Nfilt,params_fast.OutDim,params_fast.Niter);
obs_slow.init.Ytrue_full_story.val(1,1,:) = params_sim.out.simout.ECM_Vb.Data(1:tend+1)';

obs_fast.init.Nw_Nts = Nts_fast*Nw_fast;
obs_slow.init.Nw_Nts = Nts_slow*Nw_slow;


% same ground truth for obs_slow
obs_slow.init.X = obs_fast.init.X;
obs_slow.init.Ytrue_full_story = obs_fast.init.Ytrue_full_story;

% assign output
obs = {obs_fast, obs_slow};
params_out = {params_fast, params_slow};
% the whole process could be long, why not going for a nap? No worries, 
% this "sounds" like a nice way to wake up. (Uncomment)
% load handel
% sound(y,Fs)


end

