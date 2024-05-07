%% PARAMS_BATTERY
% file: params_battery.m
% author: Federico Oliva
% date: 27/05/2022
% description: this function initialises the parameters for a battery
% model (both for simulation and observation)
% INPUT: none
% OUTPUT:
% params: structure with all the necessary parameters
function params = params_battery_simulink_calce

    %%%%%%%%%%% LOAD DATA OF A BATTERY EXPERIMENT %%%%%%%%%%%

    % loading input and output signal data
    input_output_signal_filepath = '/data/CALCE/INR_18650/DST/SP2_25C_DST/11_05_2015_SP20-2_DST_50SOC.xlsx';
    [final_input_output_signal_data, params.out] = data_loader(input_output_signal_filepath);
    
    % loading ground truth SOC-OCV data
    soc_ocv_filepath = '/data/CALCE/INR_18650/Sample_1_SOC_incremental/SP1_25C_IC_OCV_12_2_2015.xlsx';
    soc_ocv = readtable(soc_ocv_filepath, 'Sheet', 'SOC_OCV');
    final_soc_ocv = soc_ocv(:, {'SOC', 'OCV'});

    % Loading input signals and parameter data
    % input_data = load('data/ECM_parameters_updated.mat');
    input_data.OCV = transpose(final_soc_ocv.OCV);
    input_data.SOC = transpose(final_soc_ocv.SOC);
    input_data.R0 = ones(size(input_data.OCV)) * 0.015;
    input_data.R1 = ones(size(input_data.OCV)) * 0.0760;
    input_data.C1 = ones(size(input_data.OCV)) * 3000;
    params.input_data = input_data;

    % time
    params.Ts = 1e0;
    params.time = final_input_output_signal_data.Step_Time_s_;

    % handle SOC not zero
    params.input_data.SOC(find(params.input_data.SOC==0)) = eps;
    
    % GET MIN AND MAX OF OCV (MEASURE)
    params.min_params = min([input_data.OCV;input_data.R0;input_data.R1;input_data.C1],[],2);
    params.max_params = max([input_data.OCV;input_data.R0;input_data.R1;input_data.C1],[],2);

    % SETUP THE EXPERIMENT - MODEL PERTURBARION
    params.deltaModel = 0*0.05;
    
    % SETUP THE EXPERIMENT - NOMINAL DATA    
    npoints = length(params.input_data.OCV);
    params.input_data.OCV_nominal = params.input_data.OCV.*(1+0.1*params.deltaModel*randn(1,npoints));
    params.input_data.R0_nominal = params.input_data.R0.*(1+params.deltaModel*randn(1,npoints));
    params.input_data.R1_nominal = params.input_data.R1.*(1+params.deltaModel*randn(1,npoints));
    params.input_data.C1_nominal = params.input_data.C1.*(1+params.deltaModel*randn(1,npoints));    
           
    % SETUP THE EXPERIMENT  - Battery Capacity (converting Ampere-hour to Ampere-second)
    params.InputAmplitude = -1;
    params.C_n_h_nominal = 2.0027*abs(params.InputAmplitude);
    params.C_n_nominal = params.C_n_h_nominal * 3600;         

    % SETUP THE EXPERIMENT - generate modular HPPC
    % define the input current module
    % params.input_current_Ts = 1;        
    % [HCCP, tspan, tspan_pos] = generate_HCCP(params.input_current_Ts,params.C_n_h_nominal);
    % params.startpos = tspan_pos(1);
    % params.stoppos = tspan_pos(end); 
    % params.input_current = HCCP;    
    % params.input_current_modular_period = params.stoppos-params.startpos;    
    % params.input_current_modular_time = 0:params.input_current_modular_period;
    % params.input_current_modular_time_dense = 0:params.input_current_Ts:params.input_current_modular_period;    
    % params.input_current_modular = interp1(params.input_current_modular_time,params.input_current(params.startpos:params.stoppos),params.input_current_modular_time_dense);
    % 
    % % slow modular HPPC - dense realization
    % params.time_slow = 1;
    % params.input_current_modular_period_slown = params.input_current_modular_period*params.time_slow;
    % params.input_current_modular_time_slown_dense = 0:params.input_current_Ts:params.input_current_modular_period_slown;
    % for i=1:length(params.input_current_modular_time_dense)
    %     params.input_current_modular_slown(params.time_slow*(i-1)+1:params.time_slow*(i-1)+params.time_slow) = params.input_current_modular(i);
    % end    

    % noise characteristics
    noise = 1;
    params.percNoise = noise*5e-2;
    params.NoisePwr = noise*5e-3;

    % temperature
    params.Temperature = 298.15;            
    
    % state dimension
    params.dim_state = 30;    
    params.dim_out = 1;
    params.dim_input = 1;
    
    % initial condition
    params.x0_simulink = zeros(params.dim_state,1);

    params.SIMstate = timeseries(zeros(params.dim_state,numel(params.time)),params.time);
    params.SIMmeasure = timeseries(zeros(params.dim_out,numel(params.time)),params.time);
    params.SIMinput = timeseries(zeros(params.dim_input,numel(params.time)),params.time);
    
end