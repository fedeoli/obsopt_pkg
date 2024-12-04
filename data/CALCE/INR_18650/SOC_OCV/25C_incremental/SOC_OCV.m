clear all;

% Load the Excel file
file_path = "/Users/tkddesai/surfdrive/Main drive/Projects/MatLab/obsopt_pkg/data/CALCE/INR_18650/SOC_OCV/25C_incremental/12_2_2015_Incremental OCV test_SP20-3.xlsx";

% Get sheet names
[~, sheet_names] = xlsfinfo(file_path);

% Exclude the 'Info' sheet
sheet_names(strcmp(sheet_names, 'Info')) = [];

% Initialize variables to store rows before step index changes to 5
rows_before_step_change_all_sheets = {};

% Loop through each sheet
for i = 1:length(sheet_names)
    % Read the data from the current sheet
    data = readtable(file_path, 'Sheet', sheet_names{i});
    
    % Extract relevant columns
    step_index = data.Step_Index;
    
    % Find the indices where the step index changes to 5 and is not already 5
    step_change_indices = find((step_index == 5 |step_index == 9) & ([0; diff(step_index)] ~= 0));
    
    % Get the row just before each step index change to 5
    rows_before_step_change = step_change_indices - 1;
    
    % Handle cases where the index might be out of bounds
    rows_before_step_change(rows_before_step_change < 1) = [];
    
    % Save the rows before step index changes to 5
    rows_before_step_change_all_sheets{i} = data(rows_before_step_change, :);
end

% Combine all rows before step index changes to 5 into one table
all_rows_before_step_change = vertcat(rows_before_step_change_all_sheets{:});

% Display the results
disp('Rows before step index changes to 5:');
disp(all_rows_before_step_change);

% Save to a new Excel file
output_file_path = 'rows_before_step_changes_3.xlsx';
writetable(all_rows_before_step_change, output_file_path);