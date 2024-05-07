function [final_input_output_signal_data, simout_struct] = data_loader(input_output_signal_filepath)

    input_output_signal_data = readtable(input_output_signal_filepath, 'Sheet','Channel_1-008');

    % changing the sign associated with current as we consider +ve sign as
    % discharge and negative sign as a charge
    input_output_signal_data.Current_A_ = -1 * input_output_signal_data.Current_A_;

    % extract rows where Step_Index turns 7, DST discharge part starts from
    % 50% SOC, also not taking all dataset into account
    startIndex = find(input_output_signal_data.Step_Index == 7, 1, 'first');
    extracted_data = input_output_signal_data(startIndex:end, :);

    % extract the columns 'Current_A_', 'Date_Time', and 'Step_Time_s_' from the extracted rows
    final_input_output_signal_data= extracted_data(:, {'Current_A_', 'Date_Time', 'Step_Time_s_', 'Voltage_V_'});
    nRows = size(final_input_output_signal_data, 1);
    % timeInSeconds = seconds(final_input_output_signal_data.Date_Time - final_input_output_signal_data.Date_Time(1));
    timeInSeconds = 0:(nRows-1);
    final_input_output_signal_data.Step_Time_s_ =timeInSeconds';

    % Create timeseries objects
    ECM_Vb = timeseries(final_input_output_signal_data.Voltage_V_, final_input_output_signal_data.Step_Time_s_, 'Name', 'ECM_Vb');
    u = timeseries(final_input_output_signal_data.Current_A_, final_input_output_signal_data.Step_Time_s_, 'Name', 'u');

    % no additional noise for now
    ECM_Vb_noise = ECM_Vb;
    
    % Capacity of the battery in Ampere-hours
    capacity = 2.0027;  % Given as 2 Ah
    
    % Calculate SOC for each time step
    soc  = 1- (extracted_data.Discharge_Capacity_Ah_ - (extracted_data.Charge_Capacity_Ah_ - extracted_data.Charge_Capacity_Ah_(1)))/capacity;
    
    % Add the SOC column to your data table
    final_input_output_signal_data.SOC = soc;
    ECM_soc = timeseries(soc, final_input_output_signal_data.Step_Time_s_, 'Name','ECM_soc');


    simout = struct('ECM_Vb_noise', ECM_Vb_noise, ...
                'ECM_Vb', ECM_Vb, ...
                'u', u, ...
                'ECM_soc', ECM_soc);
    simout_struct = struct('simout', simout,'tout', final_input_output_signal_data.Step_Time_s_);

end