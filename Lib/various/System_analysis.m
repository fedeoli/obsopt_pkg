%% script for the theoretical test

% stuff
close all
clear 
clc

% down sampling
N = 2;
plot_flag = 1;
filter = 1;
params = 1;

% sym vars
syms y [1 N] real
syms y_dot [1 N] real
syms x0 real
syms theta real
syms T real
syms V real
syms k real

% model
theta_true = 1;
x0_true = 1.1052;y
T_true = 0.1;
if params
    state = [x0, theta];
else
    state = [x0];
end

% flow
phi = x0.*exp(theta*k*T);
phi_dot = k*x0.*theta*exp(theta*k*T);
phi_fun = symfun(phi,[x0, theta, k, T]);
phi_fun_dot = symfun(phi_dot,[x0, theta, k, T]);

% check pivots around y
y_star = zeros(1,N);
for i=1:N
    y_star(i) = double(phi_fun(x0_true, theta_true, i-1,T_true));
    y_star_dot(i) = double(phi_fun_dot(x0_true, theta_true, i-1,T_true));
end

% define function
V = sym(0);
for i=1:N
   eval(['tmp = (y',num2str(i),' - phi_fun(x0,theta,',num2str(i-1),',T))^2;']);
   V = V + tmp;
   
   if filter
       eval(['tmp_dot = (y_dot',num2str(i),' - phi_fun_dot(x0,theta,',num2str(i-1),',T))^2;']);
       V = V + tmp_dot;
   end
end
if filter
    V_fun = symfun(V,[x0,theta,T,y,y_dot]);
else
    V_fun = symfun(V,[x0,theta,T,y]);
end

% compute gradient
V_grad = gradient_sym(V,state);
V_grad = simplify(V_grad);

% hessian sym
V_hess = gradient_sym(V_grad,state);
V_hess = simplify(V_hess);
if filter
    V_hess_fun = symfun(V_hess,[x0,theta,T,y,y_dot]);
else
    V_hess_fun = symfun(V_hess,[x0,theta,T,y]);
end


% function evaluation
% define evaluation grid
Tgran = [5e-2;5e-2];
x_range = [0.8*params.X, 1.2*params.X];
for i=1:params.StateDim
    x_grid(i).val = x_range(i,1):Tgran(i):x_range(i,2);
end
x0_range = x_range(1).val;
if params
    theta_range = x_range(2).val;
else
    theta_range = theta_true*ones(1,2);
end
x0_grid = x0_range(1)*x0_true:Ts:x0_range(2)*x0_true;
theta_grid = theta_range(1)*theta_true:Ts:theta_range(2)*theta_true;
[X,THETA] = meshgrid(x0_grid,theta_grid);
F = {};
F_tot = cell(size(V_hess,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot V
tmp_str = 'V_val = double(V_fun(X,THETA,T_true';
for i=1:N
   tmp_str = [tmp_str,',y_star(',num2str(i),')'];  
end
if filter
   for i=1:N
       tmp_str = [tmp_str,',y_star_dot(',num2str(i),')'];  
    end 
end
tmp_str = [tmp_str, '));'];
eval(tmp_str);
if plot_flag
    % plot
    figure;
    if length(state) > 1
        surf(X,THETA,V_val);
    else
        plot3(X,THETA,V_val);
        grid on
    end

    % stuff
    xlabel('x_0');
    ylabel('\theta');
    zlabel('V');
    title('V');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot V_hess eig
tmp_str = 'V_hess_val = V_hess_fun(X,THETA,T_true';
for i=1:N
   tmp_str = [tmp_str,',y_star(',num2str(i),')']; 
end
if filter
   for i=1:N
       tmp_str = [tmp_str,',y_star_dot(',num2str(i),')'];  
    end 
end
tmp_str = [tmp_str, ');'];
eval(tmp_str);

for i=1:length(x0_grid)
   for j=1:length(theta_grid)
       if params
           for a=1:length(state)
               for b=1:length(state)
                   tmp_V(a,b) = double(V_hess_val{a,b}(i,j));
               end
           end
       else
           tmp = V_hess_val';
           tmp_V = double(tmp(i,j));
       end
       V_hess_eig(i,j,:) = eig(tmp_V); 
   end
end

if plot_flag
    for i=1:size(V_hess,1)
        figure()
        EIG = reshape(V_hess_eig(:,:,i),[length(x0_grid),length(theta_grid)]);
        
        if length(state) > 1
            surf(X,THETA,EIG);
        else
            plot3(X,THETA,EIG);
            grid on
        end
        
        % stuff
        xlabel('x_0');
        ylabel('\theta');
        zlabel('EIG');
        title(['EIG',num2str(i)]);
        
    end
end

for i=1:size(V_hess,1)
   V_hess_eig_min(i) = min(min(V_hess_eig(:,:,i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%