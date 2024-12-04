function [final_input_output_signal_data, simout_struct] = data_loader(input_output_signal_filepath)

    % Read the specified sheet from the Excel file into a table
    input_output_signal_data = readtable(input_output_signal_filepath, 'Sheet','Channel_1-008');

    % Change the sign of the current to consider positive as discharge and negative as charge
    input_output_signal_data.Current_A_ = -1 * input_output_signal_data.Current_A_;

    % Find the starting index where Step_Index equals 7
    startIndex = find(input_output_signal_data.Step_Index == 7, 1, 'first');
    extracted_data = input_output_signal_data(startIndex:end, :);

    % Select relevant columns for further processing
    final_input_output_signal_data = extracted_data(:, {'Current_A_', 'Date_Time', 'Step_Time_s_', 'Voltage_V_'});
    nRows = size(final_input_output_signal_data, 1);

    % Generate a time vector in seconds
    timeInSeconds = 0:(nRows-1);
    final_input_output_signal_data.Step_Time_s_ = timeInSeconds';

    % Define battery parameters
    params.InputAmplitude = -1;
    params.C_n_h_nominal = 2.0 * abs(params.InputAmplitude);  % Ampere-hour
    params.C_n_nominal = params.C_n_h_nominal * 3600;        % Convert to Ampere-second

    % Calculate SOC for each time step
    soc = 0.8 - cumsum(final_input_output_signal_data.Current_A_ / (3600 * params.C_n_h_nominal));

    % Add the SOC column to the data table
    final_input_output_signal_data.SOC = soc;

    % Create a timeseries object for SOC
    ECM_soc = timeseries(soc, final_input_output_signal_data.Step_Time_s_, 'Name','ECM_soc');

    %% **Filtering Step: Exclude Data with SOC < 0.09 Based on ECM_soc.data**
    
    % Identify rows where SOC is greater than or equal to 0.09
    valid_idx = ECM_soc.Data >= 0.09;
    
    % Filter the data table to include only valid rows
    final_input_output_signal_data = final_input_output_signal_data(valid_idx, :);
    
    % Optional: Reset the time vector after filtering (if needed)
    % final_input_output_signal_data.Step_Time_s_ = 0:(height(final_input_output_signal_data)-1)';
    
    % Update the ECM_soc timeseries to reflect the filtered SOC
    ECM_soc = timeseries(final_input_output_signal_data.SOC, final_input_output_signal_data.Step_Time_s_, 'Name','ECM_soc');

    %% **Update Timeseries Objects with Filtered Data**
    
    % Create timeseries objects for Voltage and Current with the filtered data
    ECM_Vb = timeseries(final_input_output_signal_data.Voltage_V_, final_input_output_signal_data.Step_Time_s_, 'Name', 'ECM_Vb');
    u = timeseries(final_input_output_signal_data.Current_A_, final_input_output_signal_data.Step_Time_s_, 'Name', 'u');

    % For now, no additional noise is added
    ECM_Vb_noise = ECM_Vb;
    
    %% **Prepare the Output Structure for Simulink**
    
    simout = struct('ECM_Vb_noise', ECM_Vb_noise, ...
                    'ECM_Vb', ECM_Vb, ...
                    'u', u, ...
                    'ECM_soc', ECM_soc);
    
    simout_struct = struct('simout', simout, 'tout', final_input_output_signal_data.Step_Time_s_);

end
