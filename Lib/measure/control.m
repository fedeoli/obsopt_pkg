%% measure function
function u = control(x,params)

    % control law
    u = [params.K1, params.K2, params.K3]*x(1:3);
    
    % input enable
    u = params.input_enable*u;
end