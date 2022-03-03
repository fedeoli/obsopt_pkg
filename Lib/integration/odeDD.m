%% runge kutta 4 integration
function out = odeDD(ode,tspan,y0,options,varargin)

%%% general stuff %%%
solver_name = 'ode45';

% Check inputs
if nargin < 4
  options = [];
  if nargin < 3
    y0 = [];
    if nargin < 2
      tspan = [];
      if nargin < 1
        error(message('MATLAB:ode45:NotEnoughInputs'));
      end  
    end
  end
end

% Output
FcnHandlesUsed  = isa(ode,'function_handle');

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
  options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType] = ...
  odearguments_edit(FcnHandlesUsed, solver_name, ode, tspan, y0, options, varargin);


N = length(y0);             % number of state's components
M = length(tspan);         % numer of time steps
dt = tspan(2) - tspan(1);   % time step

% Matrices allocation
X = zeros(N,M);
X(:,1) = y0;
K1 = zeros(N,M);

for i = 1:M-1
    
    % State and time at ti
    x = X(:,i);
    
    % Runge Kutta 4
    K1(:,i) = feval(ode,t0,y0,odeArgs{:});
    
    % Solution at ti+1
    X(:,i+1) = K1(:,i);
    
end

% store
out.y = X;


end