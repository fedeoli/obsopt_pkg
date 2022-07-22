%% MEASURE_BATTERY
% file: measure_battery.m
% author: Federico Oliva
% date: 27/05/2022
% description: this function implements the output mapping of a battery
% model
% INPUT:
% x: state vector
% params: structure with all the necessary parameters 
% t: time instant (may be not used)
% OUTPUT:
% y: output measurement
function y = measure_battery_tushar(x,params,t)

    % compute the control - I
    params.u = params.input(t,x,params);    
    
    % get the observed components of the state vector
    y = params.Voc - x(2,:) -(params.R0 * params.u(1,:));
end