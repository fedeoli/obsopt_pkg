%% PARAMS_OSCILLATOR_VDP
% file: params_control_test.m
% author: Federico Oliva
% date: 22/06/2022
% description: this function initialises the parameters for an unstable LTI
% model
% INPUT: none
% OUTPUT:
% params: structure with all the necessary parameters
function params = params_TCV_Zaccarian(varargin)

    % case 1 = strong
    % case 2 = static annihilator
    % case 3 = dyn annihilator
    params.case = 3;
    
    % Zaccarian
    % generate system
    params.A_def = [-0.157 -0.094; ...
                -0.416 -0.45];
    params.B_def = [0.87 0.253 0.743; ...
                0.39 0.354 0.65];
    params.C_def = [0 1];
    params.D_def = [0 0 0];
    params.sys_P_def = ss(params.A_def,params.B_def,params.C_def,params.D_def);
    
    % generate dimensions
    params.n = size(params.A_def,1);
    params.m = size(params.B_def,2);
    params.p = size(params.C_def,1);
    
    % companion form
    params.A = [-0.607, 1;
                      -0.03155, 0];
    params.B = [0.39 0.354 0.65;
                      -0.3007 -0.04967 -0.207];
    params.C = [1 0];
    params.D = [0 0 0]; 
    params.sys_P = ss(params.A, params.B, params.C, params.D);       

    % generate subsystem bar
    params.q_pos = 1;    
    params.q = length(params.q_pos); 
    params.eta = [];
    
    % reference signal system (Sigma_r)
    params.dim_state_r = size(params.C,1);
    params.Ar = zeros(params.dim_state_r);
    params.Br = 0;
    params.Cr = zeros(params.dim_state_r);
    params.Dr = 1;   
    params.sys_R = ss(params.Ar,params.Br,params.Cr,params.Dr);    
    params.sys_R.InputName = 'ref';
    params.sys_R.OutputName = 'r';
    params.sys_Sum = sumblk('e = r - y',params.q);
                    
    % input stuff
    params.dim_input = params.m;    
    
    % controller system (Sigma_c)    
    params.dim_state_c = 4;
    params.Ac = [-1.57 0.5767 0.822 -0.65;...
                 -0.9 -0.501 -0.94 0.802; ...
                 0 1 -1.61 1.614; ...
                 0 0 0 0];
    params.Bcr = [0; 0; 0; 1];
    params.Bc = -params.Bcr;
    params.Cc = [1.81 -1.2 -0.46 0; ...
                 -0.62 1.47 0.89 0; ...
                 0 0 0 0];
    params.Dcr = [0; 0; 0];
    params.Dc = -params.Dcr;
    params.sys_C = ss(params.Ac,[params.Bc params.Bcr],params.Cc,[params.Dc params.Dcr]); 
    
    % only one input
    params.sys_C_err = ss(params.Ac,params.Bcr,params.Cc,params.Dcr);   
    params.sys_C_err.InputName = 'e';
    params.sys_C_err.OutputName = 'yc';
    
    %%% define PSI    
    Npoints = 5;
    poles = -logspace(0,1,Npoints);

    for i=1:Npoints-1 
        params.POLES(i,:) = linspace(poles(i),poles(i+1),3);
        params.PSI(i,:) = poly(params.POLES(i,:));
        params.PSI(i,:) = params.PSI(i,:)./params.PSI(i,end);    
    end       
%     params.tpoints = CustomStartPointSet(params.PSI(:,1:end-1));
    %%%%%%%%%%%%%%%%%%
    
    % optimizer (Sigma_op) 
    if params.case == 1
        % input allocation Strong (Zaccarian)
        params.Bperp = null([params.B; params.D]);
        params.dim_state_op = size(params.Bperp,2);
        
        % Sigma op
        params.R = [1 1 1].*eye(params.dim_input);
        params.gamma = 1e-1;
        params.A_op = -params.gamma*params.Bperp'*params.R*params.Bperp;
        params.B_op = -params.gamma*params.Bperp'*params.R;        
        params.C_op = eye(params.dim_state_op);
        params.D_op = zeros(params.dim_state_op,params.m);
        params.sys_op = ss(params.A_op,params.B_op,params.C_op,params.D_op);
        
        % Sigma an
        params.dim_state_an = 1;
        params.A_an = zeros(params.dim_state_an);
        params.B_an = zeros(params.dim_state_an,params.dim_state_op);
        params.C_an = zeros(params.m,params.dim_state_an);
        params.D_an = params.Bperp;
        params.sys_An = ss(params.A_an,params.B_an,params.C_an,params.D_an);
        
        params.NumPsi = 0;
        params.Psi = [];

    elseif params.case == 2
        % input allocation weak static (Zaccarian)  
        params.Pstar = dcgain(params.sys_P);
        params.Pperp = null(params.Pstar);        
        params.dim_state_op = size(params.Pperp,2);
        
        % Sigma op
        params.R = [1 1 1].*eye(params.dim_input);
        params.gamma = 1e-1;
        params.A_op = -params.gamma*params.Pperp'*params.R*params.Pperp;
        params.B_op = -params.gamma*params.Pperp'*params.R;        
        params.C_op = eye(params.dim_state_op);
        params.D_op = zeros(params.dim_state_op,params.m);
        params.sys_op = ss(params.A_op,params.B_op,params.C_op,params.D_op);
        
        % Sigma an
        params.dim_state_an = 1;
        params.A_an = zeros(params.dim_state_an);
        params.B_an = zeros(params.dim_state_an,params.dim_state_op);
        params.C_an = zeros(params.m,params.dim_state_an);
        params.D_an = params.Pperp;
        params.sys_An = ss(params.A_an,params.B_an,params.C_an,params.D_an);
        
        params.NumPsi = 0;
        params.Psi = [];
        
    else
        % input allocation weak dynamic (Zaccarian)  
        % dynamic annihilator    
        [params.num_An,params.N, params.N_An] = annihilator(params.sys_P, params.n + 1); 
        W_An = tf(params.num_An, params.PSI(end,:));
        % Compute the annihilator state-space model
        params.sys_An = ss(W_An);
        params.A_an = params.sys_An.A;
        params.B_an = params.sys_An.B;
        params.C_an = params.sys_An.C;
        params.D_an = params.sys_An.D;
        params.dim_state_an = size(params.A_an,1);
        params.Psi = params.PSI(1,:)';
%         params.Psi = [2.7118e+03   1.0795e+06   1.0781e+08   1.0000e+00]';
        params.NumPsi = length(params.Psi);                            
                
        params.Anstar = dcgain(params.sys_An);        
        params.dim_state_op = size(params.Anstar,2);
        
        % Sigma op
        %%% first realisation (non minimal)
        params.R = [1 1 1].*eye(params.dim_input);        
        params.A_op_def = -params.Anstar'*params.R*params.Anstar;
        params.B_op_def = -params.Anstar'*params.R;        
        params.C_op_def = eye(params.dim_state_op);
        params.D_op_def = zeros(params.dim_state_op,params.m);
%         params.sys_op_def = minreal(ss(params.A_op_def,params.B_op_def,params.C_op_def,params.D_op_def),[],false);
        params.sys_op_def = ss(params.A_op_def,params.B_op_def,params.C_op_def,params.D_op_def);
        %%% minimal realisation (update with GAMMA)
        params.dim_state_op = size(params.sys_op_def.A,2);
        
        %%% define GAMMA    
        Npoints = 5;
        decades = logspace(-2,1,Npoints);

        for i=1:Npoints 
            params.GAMMA(i,:) = repmat(decades(i),1,params.dim_state_op);            
        end      
        % define tpoints        
        params.tpoints = CustomStartPointSet(params.GAMMA);
        %%%%%%%%%%%%%%%%%%
        
%         params.gamma = params.GAMMA(1,:)';
        params.gamma = 1.0*ones(params.dim_state_op,1);       
%         params.gamma = 1e0*[1.9429e+02   4.1329e+02   1.4146e+00   1.1411e+00   8.5913e+01]';
        params.Gamma = diag(params.gamma);
        params.NumGamma = params.dim_state_op;
        params.A_op = params.Gamma*params.sys_op_def.A;
        params.B_op = params.Gamma*params.sys_op_def.B;
        params.C_op = params.sys_op_def.C;
        params.D_op = params.sys_op_def.D;
        params.sys_op = ss(params.A_op,params.B_op,params.C_op,params.D_op);                    
                
        % get steady state input 
        A_inf = [params.Ac, -params.Bc * params.C;
                params.B*params.Cc, params.A];
        B_inf = [params.Bc;
                 zeros(params.n, params.p)];
        C_inf = [params.Cc, zeros(params.m, params.n)];
        y_c_inf = -C_inf * pinv(A_inf) * B_inf;
        params.u_inf = (eye(params.m) - params.Anstar * pinv(params.Anstar' * params.R * params.Anstar) * params.Anstar' * params.R) * y_c_inf;
        
    end    
    
    params.sys_An.InputName = 'v';
    params.sys_An.OutputName = 'ya';
    params.sys_op.InputName = 'yc';
    params.sys_op.OutputName = 'v';
    
                            
    % state dimension
    params.dim_state = params.dim_state_r + params.dim_state_c + params.dim_state_op + params.dim_state_an + params.n + params.NumPsi + params.NumGamma;        
    
    % initial condition
    % [x, xc, xop, xan, psi]
    params.X(1).val(:,1) = [zeros(params.dim_state_r,1); zeros(params.dim_state_c,1); zeros(params.dim_state_op,1); zeros(params.dim_state_an,1); 0*ones(params.n,1); params.Psi; params.gamma];
    
    % position in the state vector of the estimated parameters
    params.estimated_params = [params.dim_state_r + params.dim_state_c + params.dim_state_op + params.dim_state_an + params.n + 1:params.dim_state];
    
    % which vars am I optimising
    %%% PSI opt %%%
    %params.opt_vars = [params.dim_state_r + params.dim_state_c + params.dim_state_op + params.dim_state_an + params.n + 1:params.dim_state-params.NumGamma - 1];
    %%% GAMMA opt %%%
    %params.opt_vars = [params.dim_state-params.NumGamma+1:params.dim_state];
    %%% ALL OPT %%%
    params.PsiPos = params.dim_state_r + params.dim_state_c + params.dim_state_op + params.dim_state_an + params.n + 1:params.dim_state-params.NumGamma;
    params.GammaPos = params.dim_state-params.NumGamma+1:params.dim_state;
    params.opt_vars = [params.PsiPos(1:end-1) params.GammaPos];
    
    % set the not optimised vars
    tmp = 1:length(params.X(1).val(:,1));
    tmp_idx = tmp;
    for i=1:length(params.opt_vars)
        tmp_idx = intersect(tmp_idx,find(tmp~=params.opt_vars(i)));
    end
    params.nonopt_vars = tmp_idx;
    
    % out vars
    params.OutDim = 1;
    params.OutDim_compare = [1];
    
    % plot vars (used to plot the state estimation. When the parameters are
    % too many, consider to use only the true state components)
    params.plot_vars = 1:(params.n + params.dim_state_c + params.dim_state_r + params.dim_state_op + params.dim_state_an);
    params.plot_params = (params.n + params.dim_state_c + params.dim_state_r + params.dim_state_op + params.dim_state_an + 1):params.dim_state;   
    
    % number of reference trajectories (under development)
    params.Ntraj = 2;
    params.traj = 1;
    params.optimising = 0;
    params.Ru = 0.25e-1;
    
    % perturbed models
    params.sys_pert(1).A = params.A;
    params.sys_pert(1).B = params.B;
    params.sys_pert(1).C = params.C;
    params.sys_pert(1).D = params.D;
    params.sys_pert(1).sys_P = params.sys_P;
    
    % pert perc
    params.pert_perc = 0.1;
    for i=2:params.Ntraj
        %%% no change in Pstar
%         params.sys_pert(i).A = [[params.A(1:end-1,1).*(1+params.pert_perc*randn(params.n-1,1)); params.A(end,1)] params.A(:,2:end)];
%         params.sys_pert(i).B = [params.B(1:end-1,:).*(1+params.pert_perc*randn(params.n-1,params.m)); params.B(end, :)];
        %%% chage Pstar
        params.sys_pert(i).A = [params.A(1:end,1).*(1+params.pert_perc*rand(params.n,1)) params.A(:,2:end)];
        params.sys_pert(i).B = params.B.*(1+params.pert_perc*rand(params.n,params.m));
        params.sys_pert(i).C = params.C;
        params.sys_pert(i).D = params.D;
        params.sys_pert(i).sys_P = ss(params.sys_pert(i).A,params.sys_pert(i).B,params.sys_pert(i).C,params.sys_pert(i).D);
        
        params.X(i).val(:,1) = params.X(i-1).val(:,1);
    end
    
    %%% closed loop system %%%
    for i=1:params.Ntraj
        %%% no allocation closed loop %%%
        % plant and controller        
        params.sys_pert(i).sys_P.InputName = 'u';
        params.sys_pert(i).sys_P.OutputName = 'y';        
        params.sys_SumAll = sumblk('u = yc',params.m);        
        params.sys_pert(i).sys_CLu = connect(params.sys_Sum,params.sys_C_err,params.sys_SumAll,params.sys_pert(i).sys_P,'r','u');
        params.sys_pert(i).sys_CL = connect(params.sys_Sum,params.sys_C_err,params.sys_SumAll,params.sys_pert(i).sys_P,'r','y');         
        
        params.sys_SumAll = sumblk('u = yc + ya',params.m);  
        params.sys_pert(i).sys_CL_Allu = connect(params.sys_Sum,params.sys_C_err,params.sys_op,params.sys_An,params.sys_SumAll,params.sys_pert(i).sys_P,'r','u');
        params.sys_pert(i).sys_CL_All = connect(params.sys_Sum,params.sys_C_err,params.sys_op,params.sys_An,params.sys_SumAll,params.sys_pert(i).sys_P,'r','y');
        
    end
        
end