%% SIMULATION_GENERAL_V3
% file: simulation_general_v3.m
% author: Federico Oliva
% date: 10/01/2022
% description: function to setup and use the MHE out.observer on general model
% INPUT: none
% OUTPUT: params,out.obs
% plot results for control design
function plot_general_calce(obj,varargin)
    
    set(0,'DefaultFigureWindowStyle','docked');            
    
    fontsize = 20;
    fig_count = 0;
    
    %%%% plot state estimation %%%
    fig_count = fig_count+1;
    figure(fig_count)            
    sgtitle('State estimation')
    ylabels = {'$Z$', '$V1$', '$OCV$'};
    for i=1:length(obj.setup.plot_vars)
        subplot(length(obj.setup.plot_vars),1,i);
        hold on
        grid on
        box on
        
        for traj=1:obj.setup.Ntraj
            if 1 || strcmp(obj.setup.DataType,'simulated')
                plot(obj.setup.time,obj.init.X(traj).val(obj.setup.plot_vars(i),:),'b--','LineWidth',2);
            end
            plot(obj.setup.time,obj.init.X_est_runtime(traj).val(obj.setup.plot_vars(i),:),'r--','LineWidth',2);                                                  
        end    

        % labels
        set(gca,'fontsize', fontsize) 
        ylabel(ylabels{i},'Interpreter','latex');
    end
    % legend
    if strcat(obj.setup.DataType,'simulated')
        legend('True','Est')
    else
        legend('Stored','Est','Runtime')
    end    
    xlabel('time [s]')
    
    
    %%%% plot parameters estimation %%%
    if ~isempty(obj.setup.plot_params)
    
        % Define groups of parameters
        groups = {
            [4, 5, 6]... 
            [7, 11, 15, 19]... 
            [8, 12, 16, 20]... 
            [9, 13, 17, 21]... 
            [10, 14, 18, 22]
        };
        ylabels = {
                   ['$R0$'; '$R1$'; '$C1$'], ...
                   ['$\alpha_{0,OCV}$'; '$\alpha_{1,OCV}$'; '$\alpha_{2,OCV}$'; '$\alpha_{3,OCV}$'], ...
                   ['$\alpha_{0,R0}$'; '$\alpha_{1,R0}$'; '$\alpha_{2,R0}$'; '$\alpha_{3,R0}$'], ...
                   ['$\alpha_{0,R1}$'; '$\alpha_{1,R1}$'; '$\alpha_{2,R1}$'; '$\alpha_{3,R1}$'], ...
                   ['$\alpha_{0,C1}$'; '$\alpha_{1,C1}$'; '$\alpha_{2,C1}$'; '$\alpha_{3,C1}$']
                   };        

        % Iterate over each group
        for g = 1:length(groups)
            fig_count = fig_count + 1;
            figure(fig_count)
            sgtitle('Parameters estimation')
    
            % Plot each parameter in the group
            for i = 1:length(groups{g})
                subplot(length(groups{g}), 1, i);
                hold on
                grid on
                box on
    
                for traj = 1:obj.setup.Ntraj
                    if 1 || strcmp(obj.setup.DataType, 'simulated')
                        plot(obj.setup.time, obj.init.X(traj).val(groups{g}(i), :), 'b--', 'LineWidth', 2);
                    end
                    plot(obj.setup.time, obj.init.X_est_runtime(traj).val(groups{g}(i), :), 'g--', 'LineWidth', 2);
                end
    
                % labels
                set(gca, 'fontsize', fontsize)
                ylabel(ylabels{g}(i,:), 'Interpreter', 'latex');
                end
    
            if strcmp(obj.setup.DataType, 'simulated')
                legend('True', 'Est')
            else
                legend('Stored', 'Est', 'Runtime')
            end
            xlabel(['time [s]'])
        end
    end            
    
    %%%% plot windowed data %%%%            
    fig_count = fig_count+1;
    figure(fig_count)
    subplot(2,1,1)
    grid on
    sgtitle('Sampled measured')
    ax = zeros(1,3);
    for k=1:obj.setup.dim_out
        
        % number fo subplots depending on the output dimension
        n_subplot = obj.setup.dim_out+1;
        
        % indicize axes
        ax_index = k;
        ax(ax_index)=subplot(n_subplot,1,ax_index);
        
        % hold on on plots
        hold on
        
        % dpwn sampling instants
        WindowTime = obj.setup.time(obj.init.temp_time);
        
        for traj=1:obj.setup.Ntraj
            % plot true values
            y_meas = reshape(obj.init.Y_full_story(traj).val(1,k,:),size(obj.setup.time));
            plot(obj.setup.time,y_meas,'m:','LineWidth',2)

            % plot target values    
            try
                data = reshape(obj.init.target_story(traj).val(1,k,obj.init.temp_time),1,length(WindowTime));
                plot(WindowTime,data,'bo','MarkerSize',5);
            catch 
                disp('CHECK T_END OR AYELS CONDITION - LOOKS LIKE NO OPTIMISATION HAS BEEN RUN')
            end            
        end
        set(gca,'fontsize', fontsize)
        ylabel(strcat('y_',num2str(k)));        
    end
    xlabel('simulation time [s]');
    legend('meas','sampled')
    linkaxes(ax(1:n_subplot-1),'x');
    
    %%% plot adaptive sampling            
    ax(n_subplot) = subplot(n_subplot,1,n_subplot);
    % frequency constraint
    % y_meas = squeeze(obj.init.Y_full_story.val(1,1,:));  
    % [WT,F] = cwt(y_meas,obj.init.wvname,1/obj.setup.Ts,'VoicesPerOctave',obj.init.Nv,'FrequencyLimits',obj.init.FLIMITS);    
    % heatmap(obj.setup.time,F,real(WT))
    % grid off
    % colormap jet
    % 
    % %%% single cwt
    % fig_count = fig_count+1;
    % figure(fig_count)
    % cwt(y_meas,obj.init.wvname,1/obj.setup.Ts,'VoicesPerOctave',obj.init.Nv,'FrequencyLimits',obj.init.FLIMITS);
            
    
end