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
    % '/data/CALCE/INR_18650/BJDST/SP2_25C_BJDST/11_12_2015_SP20-2_BJDST_80SOC.xlsx';
    % '/data/CALCE/INR_18650/US06/SP2_25C_US06/11_11_2015_SP20-2_US06_80SOC.xlsx';
    % '/data/CALCE/INR_18650/FUDS/SP2_25C_FUDS/11_06_2015_SP20-2_FUDS_80SOC.xlsx';
    % '/data/CALCE/INR_18650/DST/SP2_25C_DST/11_05_2015_SP20-2_DST_80SOC.xlsx'
    input_output_signal_filepath = '/data/CALCE/INR_18650/DST/SP2_25C_DST/11_05_2015_SP20-2_DST_80SOC.xlsx';
    [final_input_output_signal_data, params.out] = data_loader(input_output_signal_filepath);
    
    % loading ground truth SOC-OCV data
    soc_ocv_filepath = '/data/CALCE/INR_18650/Sample_1_SOC_incremental/SP1_25C_IC_OCV_12_2_2015.xlsx';
    soc_ocv = readtable(soc_ocv_filepath, 'Sheet', 'SOC_OCV');
    final_soc_ocv = soc_ocv(:, {'SOC', 'OCV'});

    % Loading input signals and parameter data
    % input_data = load('data/ECM_parameters_updated.mat');
    input_data.OCV = transpose(final_soc_ocv.OCV);
    input_data.SOC = transpose(final_soc_ocv.SOC);

    % initial parameters
    input_data.R0 = [0.02, 0.02166, 0.02166, 0.02166, 0.02166, 0.02166, 0.02333, 0.025, 0.02667, 0.02833, 0.03];
    input_data.R1 =  [0.004, 0.007, 0.007, 0.007, 0.007, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012];
    input_data.C1 = [1250, 2500, 2500, 2500, 2500, 2500, 3000, 3500, 4000, 4500, 5000];
    
    % for consistency
    params.input_data = input_data;

    % time
    %                                    
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
             
    % noise characteristics
    noise = 0;
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