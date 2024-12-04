%% SIMULATION_REALDATA_BATTERY_SIMULINK
% file: simulation_general_v3.m
% author: 
% date: 21/07/2024
% description: function to setup and use the EKF on battery from
% real data of CALCE dataset 
% INPUT: none
% OUTPUT: ekf_results (struct containing y, z, and time data)

function ekf_results = simulation_realdata_battery_EKF

    clc;
    close all;
    clear all;

    rng('default');

    % Loading the dataset
    params_sim = params_battery_simulink_calce;
    params = params_battery_calce(params_sim);
    

    % Covariances and intial SOC
    % params.initial_soc = 0.4;
    % P0 = diag([1E-0, 1E-2, 1E-6 *ones(1, 12)]); 
    % Q = diag([1E-4, 1E-3, 1E-6, 1E-6, 1E-6 *ones(1, 2), 1E-6, 1E-6, 1E-6 *ones(1, 2), 1E-6, 1E-6, 1E-6 *ones(1, 2)]);
    % R = 1E-6;

    params.initial_soc = 0.4;
    P0 = diag([1E-2, 1E-2, 1E-6 *ones(1, 12)]); 
    Q = diag([1E-4, 1E-3, 1E-6*ones(1, 12)]);
    R = 1E-6;
    
    % intermdediate
    % P0 = diag([1E-2, 1E-2, 1E-3*ones(1, 12)]); 
    % Q = diag([1E-4, 1E-4, 1E-6*ones(1, 12)]);
 
    % Coefficients
    Voc_coefficient = [params.alpha_Voc params.beta_Voc params.gamma_Voc params.delta_Voc params.eps_Voc params.xi_Voc];
    R0_coefficient = [params.alpha_R0 params.beta_R0 params.gamma_R0 params.delta_R0];
    R1_Coefficient = [params.alpha_R1 params.beta_R1 params.gamma_R1 params.delta_R1];
    C1_Coefficient = [params.alpha_C1 params.beta_C1 params.gamma_C1 params.delta_C1];


    % R0_coefficient = [0.0863 -0.1015 0.1998 -0.1199];
    % R1_Coefficient = [-0.0020 0.2940 -0.9014 0.7453];
    % C1_Coefficient = [828.894553983524 9.244656743141175e+03 -3.611686504893436e+04 3.229387762053796e+04];
    % Battery capacity (converting Ampere-hour to Ampere-second)
    capacity = params.C_n_h * 3600;

    % Sampling time
    Ts = params.Ts;

    % Input, output and SOC ground truth
    u = params.u_sim;
    y_true = params.y_true_sim;
    z_true = params.soc_sim;
    datatime = params.time;

    % Initial states
    initial_Voc = polyval(flip(Voc_coefficient), params.initial_soc);
    initial_R0 = polyval(flip(R0_coefficient), params.initial_soc);
    initial_V1 = initial_Voc - u(1)*initial_R0 - y_true(1);
    x_initial = [params.initial_soc; initial_V1; R0_coefficient'; R1_Coefficient'; C1_Coefficient'];

    % Initialize the estimated state and prediction
    x_estimated = zeros(length(x_initial), length(datatime));
    comp_time = zeros(1, length(datatime)); 
    y_predict = zeros(1, length(datatime));

    tic; 

    % EKF first step - correction based on initial estimate
    [~, H] = jacobian(x_initial, u(1), Voc_coefficient, Ts);
    K = P0 * H' / (H * P0 * H' + R);
    y_predict(1) = battery_measurement(x_initial, u(1), Voc_coefficient);
    x_estimated(:, 1) = x_initial + K * (y_true(1) - y_predict(1));
    P = (eye(length(x_initial)) - K * H) * P0;
    P = (P + P') / 2; % for symmetry

    comp_time(1) = toc; % Time taken for the first step

    for n = 2:length(datatime)
        tic;

        % Prediction
        [F, ~] = jacobian(x_estimated(:, n-1), u(n-1), Voc_coefficient, Ts);
        P = F * P * F' + Q;
        P = (P + P') / 2;
        x_predict = battery_model(x_estimated(:, n-1), u(n-1), Ts, capacity);

        % Correction
        [~, H] = jacobian(x_predict, u(n), Voc_coefficient, Ts);
        K = P * H' / (H * P * H' + R);
        y_predict(n) = battery_measurement(x_predict, u(n), Voc_coefficient);
        x_estimated(:, n) = x_predict + K * (y_true(n) - y_predict(n));
        P = (eye(length(x_initial)) - K * H) * P;
        P = (P + P') / 2;

        comp_time(n) = toc; % Time taken for each step
    end

    z_estimated = x_estimated(1, :);
    soc_rmse = sqrt(mean((z_estimated - z_true).^2));
    voltage_rmse = sqrt(mean((y_predict - y_true).^2));
    avg_comp_time = mean(comp_time); % Average computation time
    max_comp_time = max(comp_time); % Maximum computation time

    % Display RMSE and computation times
    fprintf('SOC RMSE: %.6e\n', soc_rmse);
    fprintf('Voltage RMSE: %.6e\n', voltage_rmse);
    fprintf('Mean Computational Time: %.6e seconds\n', avg_comp_time);
    fprintf('Max Computational Time: %.6e seconds\n', max_comp_time);

    % Assign results to struct
    ekf_results = struct();

    % SOC Data
    ekf_results.SOC = struct();
    ekf_results.SOC.true = z_true;
    ekf_results.SOC.estimated = z_estimated;
    ekf_results.SOC.RMSE = soc_rmse;

    % Voltage Data
    ekf_results.Voltage = struct();
    ekf_results.Voltage.true = y_true;
    ekf_results.Voltage.predicted = y_predict;
    ekf_results.Voltage.RMSE = voltage_rmse;

    % Time Data
    ekf_results.Time = datatime;

    % Computational Time Data
    ekf_results.ComputationTime = struct();
    ekf_results.ComputationTime.average = avg_comp_time;
    ekf_results.ComputationTime.maximum = max_comp_time;

    % Plot results
    plot_results(ekf_results);

end

%% Jacobian
function [F, H] = jacobian(x, u, a_OCV, Ts)
    z = x(1);
    V1 = x(2);

    R0_coeff = x(3:6)';
    R1_coeff = x(7:10)';
    C1_coeff = x(11:14)';

    soc_poly = [1; z; z^2; z^3];
    soc_deriv = [0; 1; 2*z; 3*z^2];
    soc_ocv_poly = [1; z; z^2; z^3; z^4; z^5];
    soc_ocv_deriv = [0; 1; 2*z; 3*z^2; 4*z^3; 5*z^4];

    R0 = dot(R0_coeff, soc_poly);
    R1 = dot(R1_coeff, soc_poly);
    C1 = dot(C1_coeff, soc_poly);
    tau = R1 * C1;

    R0_deriv = R0_coeff * soc_deriv; 
    R1_deriv = R1_coeff * soc_deriv; 
    C1_deriv = C1_coeff * soc_deriv; 
    tau_deriv = R1_deriv * C1 + R1 * C1_deriv; 
    OCV_deriv = a_OCV * soc_ocv_deriv; 


    % Jacobian for the state transition matrix
    F = zeros(14, 14);

    F(1, 1) = 1;

    F(2, 1) = V1 * exp(-Ts/tau) * (Ts/tau^2) * tau_deriv ...
              + u * R1_deriv * (1 - exp(-Ts/tau)) ...
              - u * R1 * exp(-Ts/tau) * (Ts/tau^2) * tau_deriv;
    F(2, 2) = exp(-Ts/tau);

    F(2, 3) = 0;
    F(2, 4) = 0;
    F(2, 5) = 0;
    F(2, 6) = 0;

    F(2, 7) = V1 * exp(-Ts/tau) * (Ts/tau^2) * C1 ...
              + u * (1-exp(-Ts/tau)) ...
              - u * R1 * exp(-Ts/tau) * (Ts/tau^2) * C1;
    F(2, 8) = V1 * exp(-Ts/tau) * (Ts/tau^2) * C1 * z ...
               + u * (1-exp(-Ts/tau)) * z ...
               - u * R1 * exp(-Ts/tau) * (Ts/tau^2) * C1 * z;
    F(2, 9) = V1 * exp(-Ts/tau) * (Ts/tau^2) * C1 * z^2 ...
               + u * (1-exp(-Ts/tau)) * z^2 ...
               - u * R1 * exp(-Ts/tau) * (Ts/tau^2) * C1 * z^2;
    F(2, 10) = V1 * exp(-Ts/tau) * (Ts/tau^2) * C1 * z^3 ...
               + u * (1-exp(-Ts/tau)) * z^3 ...
               - u * R1 * exp(-Ts/tau) * (Ts/tau^2) * C1 * z^3;

    F(2, 11) = (V1 - u * R1) * (exp(-Ts/tau) * Ts * R1/tau^2);
    F(2, 12) = (V1 - u * R1) * (exp(-Ts/tau) * Ts * R1 * z/tau^2);
    F(2, 13) = (V1 - u * R1) * (exp(-Ts/tau) * Ts * R1 * z^2/tau^2);
    F(2, 14) = (V1 - u * R1) * (exp(-Ts/tau) * Ts * R1 * z^3/tau^2);

    F(3:14, 3:14) = eye(12);

    % Jacobian for measurement matrix
    H = zeros(1, 14);

    H(1, 1) = OCV_deriv - R0_deriv * u;
    H(1, 2) = -1;
    H(1, 3) = -u;
    H(1, 4) = - u * z;
    H(1, 5) = - u * z^2;
    H(1, 6) = - u * z^3;

end

%% Battery Measurement Model
function y_true = battery_measurement(x, u, a_OCV)
    z = x(1);
    V1 = x(2);

    R0_coeff = x(3:6);

    OCV = dot(a_OCV, [1; z; z^2; z^3; z^4; z^5]);
    R0 = dot(R0_coeff, [1; z; z^2; z^3]);

    y_true = OCV - u * R0 - V1;
end

%% Battery States Model
function x_predict = battery_model(x, u, Ts, capacity)
    z = x(1);
    V1 = x(2);

    R1_coeff = x(7:10);
    C1_coeff = x(11:14);

    R1 = dot(R1_coeff, [1; z; z^2; z^3]);
    C1 = dot(C1_coeff, [1; z; z^2; z^3]);
    tau = R1 * C1;

    z = z - (u * Ts) / capacity;
    V1 = V1 * exp(-Ts / tau) + R1 * u * (1 - exp(-Ts / tau));

    x_predict = [z; V1; x(3:end)];
end

%% Plot Results
function plot_results(ekf_results)

    datatime = ekf_results.Time;
    z_true = ekf_results.SOC.true;
    z_estimated = ekf_results.SOC.estimated;
    y_true = ekf_results.Voltage.true;
    y_predict = ekf_results.Voltage.predicted;

    set(0,'DefaultFigureWindowStyle','docked');
    fontsize = 20;
    colors = [
        31, 119, 180;  % Blue
        255, 127, 14;  % Orange
        44, 160, 44;   % Green
        214, 39, 40;   % Red
        148, 103, 189; % Purple
    ] / 255;

    figure(1);
    hold on;
    grid on;
    box on;
    plot(datatime, y_true, 'LineStyle', '--', 'LineWidth', 2, 'Color', colors(1, :));
    plot(datatime, y_predict, 'LineStyle', '--', 'LineWidth', 2, 'Color', colors(2, :));
    set(gca, 'fontsize', fontsize);
    ylabel('Voltage [V]');
    xlabel('time [s]');
    legend('True Voltage', 'Predicted Voltage');

    figure(2);
    hold on;
    grid on;
    box on;
    plot(datatime, z_true, 'LineStyle', '--', 'LineWidth', 2, 'Color', colors(1, :));
    plot(datatime, z_estimated, 'LineStyle', '--', 'LineWidth', 2, 'Color', colors(2, :));
    set(gca, 'fontsize', fontsize);
    ylabel('SOC');
    xlabel('Time [s]');
    legend('True SOC', 'Estimated SOC');

end