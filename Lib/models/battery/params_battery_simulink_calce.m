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

    % input_data.R0 = [0.015,0.014167,0.013333,0.0125,0.011667,0.010833,0.010833,0.010833,0.010833,0.010833,0.01];
    % input_data.R1 = [0.024,0.022,0.02,0.018,0.016,0.014,0.014,0.014,0.014,0.014,0.008];
    % input_data.C1 = [10000,9000,8000,7000,6000,5000,5000,5000,5000,5000,2500];
    
    % input_data.R0 = [0.089000,0.076597,0.072164,0.071013,0.070835,0.070937,0.071477,0.072694,0.074148,0.073950,0.068000];
    % input_data.R1  = [0.002700,0.023089,0.016769,0.012580,0.017430,0.025575,0.027914,0.021274,0.017695,0.053725,0.199700];
    % input_data.C1 = [1877.260,1257.901,933.547,914.623,1022.945,1041.227,862.591,640.071,936.119,2872.117,8277.880];

    % main working initial parameters
    input_data.R0 = [0.01, 0.010833, 0.010833, 0.010833, 0.010833, 0.010833, 0.011667, 0.0125, 0.013333, 0.014167, 0.015]*2;
    input_data.R0 = [0.01, 0.010833, 0.010833, 0.010833, 0.010833, 0.010833, 0.011667, 0.0125, 0.013333, 0.014167, 0.015]*2;
    input_data.R1 = [0.008, 0.014, 0.014, 0.014, 0.014, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024]/2;
    input_data.R1 = [0.008, 0.014, 0.014, 0.014, 0.014, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024]/10;
    input_data.C1 = [2500, 5000, 5000, 5000, 5000, 5000, 6000, 7000, 8000, 9000, 10000]/2;

    % input_data.R0 = [0.089000,0.076597,0.072164,0.071013,0.070835,0.070937,0.071477,0.072694,0.074148,0.073950,0.068000];
    % input_data.R1 = [0.002700,0.023089,0.016769,0.012580,0.017430,0.025575,0.027914,0.021274,0.017695,0.053725,0.199700];
    % input_data.C1 = [1877.260,1257.901,933.547,914.623,1022.945,1041.227,862.591,640.071,936.119,2872.117,8277.880];
    % input_data.R0 = [0.0893053439242714,	0.0769243599952580,	0.0724835305634742,	0.0712636268626486,	0.0709228555309545,	0.0707309570101154,	0.0708033039445091,	0.0713349995802734,	0.0718349761644106,	0.0703600933438924,	0.0627492365647655];
    % input_data.R1 = [0.0049896090704552,	0.0231480479029243,	0.0169188228288767,	0.0128448105341071,	0.0178244029391147,	0.0261108491808152,	0.0286022457576631,	0.0221315266747742,	0.0187564535890487,	0.0550496059542916,	0.201388371166333];
    % input_data.C1 = [1877.25503452128,	1257.89593651177,	933.542371954046,	914.618441142885,	1022.94016538849,	1041.22249223350,	862.586300670110,	640.065406357061,	936.113566836751,	2872.11148675225,	8277.87382306438];
    %from other
    % input_data.R0 = [0.0866, 0.0781, 0.0731, 0.0706, 0.0700, 0.0705, 0.0715, 0.0721, 0.0717, 0.0694, 0.0646];
    % input_data.R1 = [-0.0028, 0.0194, 0.0273, 0.0256, 0.0191, 0.0123, 0.0101, 0.0171, 0.0381, 0.0776, 0.1404];
    % input_data.C1 = [868.393224670002, 1457.24804139726, 1511.20562111728, 1227.68298095542, 804.097138037040, 437.865109487500, 326.403912432159, 667.130563996378, 1657.46208130552, 3494.81548148494, 6376.60778166000];

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