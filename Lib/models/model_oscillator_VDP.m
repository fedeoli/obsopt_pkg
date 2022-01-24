%% mock up model
function x_dot = model_oscillator_VDP(t,x,params)

    x_dot = zeros(length(x),1);
    
    x_dot(1) = x(2);
    x_dot(2) = -x(1) + params.mu*(1-x(1)^2)*x(2);
end