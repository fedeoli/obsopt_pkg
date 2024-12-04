
% Define file paths for the four datasets
filepaths = {
    '/data/CALCE/INR_18650/BJDST/SP2_25C_BJDST/11_12_2015_SP20-2_BJDST_80SOC.xlsx';
    '/data/CALCE/INR_18650/US06/SP2_25C_US06/11_11_2015_SP20-2_US06_80SOC.xlsx';
    '/data/CALCE/INR_18650/FUDS/SP2_25C_FUDS/11_06_2015_SP20-2_FUDS_80SOC.xlsx';
    '/data/CALCE/INR_18650/DST/SP2_25C_DST/11_05_2015_SP20-2_DST_80SOC.xlsx'
};

% Define sheet names for each dataset (if different, update accordingly)
sheet_name = 'Channel_1-008';

% Define labels for the subplots
labels = {'BJDST', 'US06', 'FUDS', 'DST'};

figure;

for i = 1:length(filepaths)
    % Read data from each file
    input_output_signal_data = readtable(filepaths{i}, 'Sheet', sheet_name);

    % Find the starting index of step 7
    startIndex = find(input_output_signal_data.Step_Index == 7, 1, 'first');
    extracted_data = input_output_signal_data(startIndex:end, :);

    % Extract Date_Time and Current_A_
    Date_Time = extracted_data.Date_Time;
    Current_A_ = extracted_data.Current_A_;
    
    % Convert Date_Time to seconds from the start
    Date_Time = seconds(Date_Time - Date_Time(1));
    
    % Plotting
    subplot(2, 2, i);
    plot(Date_Time, Current_A_, 'LineWidth', 1); % Plot with thicker line
    xlabel('time [s]');
    ylabel('current [A]');
    % title(labels{i});
    grid on;
    xlim([0, max(Date_Time)]); % Ensuring x-axis starts from 0
end
