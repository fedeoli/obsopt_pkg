%% SIMULATION_GENERAL_V3
% file: simulation_general_v3.m
% author: Federico Oliva
% date: 10/01/2022
% description: function to setup and use the MHE out.observer on general model
% INPUT: none
% OUTPUT: params,out.obs
% plot results for control design
function plot_general_calce_tushar(obj,varargin)
    
    set(0,'DefaultFigureWindowStyle','docked');            
    
    fontsize = 20;
    fig_count = 0;
    category20 = [
    31, 119, 180; % Color 1: Blue
    255, 127, 14; % Color 2: Orange
    44, 160, 44; % Color 3: Green
    214, 39, 40; % Color 4: Red
    148, 103, 189; % Color 5: Purple
    140, 86, 75; % Color 6: Brown
    227, 119, 194; % Color 7: Pink
    127, 127, 127; % Color 8: Gray
    188, 189, 34; % Color 9: Olive
    23, 190, 207; % Color 10: Cyan
    174, 199, 232; % Color 11: Light Blue
    255, 187, 120; % Color 12: Light Orange
    152, 223, 138; % Color 13: Light Green
    255, 152, 150; % Color 14: Light Red
    197, 176, 213; % Color 15: Lavender
    196, 156, 148; % Color 16: Tan
    247, 182, 210; % Color 17: Light Pink
    199, 199, 199; % Color 18: Light Gray
    219, 219, 141; % Color 19: Light Olive
    158, 218, 229  % Color 20: Light Cyan
] / 255; 
    
    %%%% plot state estimation %%%
    for i=1:length(obj.setup.plot_vars)
        fig_count = fig_count+1;
        figure(fig_count)

        % First subplot: True and estimated values
        subplot(2, 1, 1);
        hold on;
        grid on;
        box on;
        true_values = [];
        estimated_values = [];
        
        for traj = 1:obj.setup.Ntraj
            if 1 || strcmp(obj.setup.DataType,'simulated')
                plot(obj.setup.time, obj.init.X(traj).val(obj.setup.plot_vars(i), :), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'blue');
                true_values = [true_values; obj.init.X(traj).val(obj.setup.plot_vars(i), :)];
            end
            plot(obj.setup.time, obj.init.X_est_runtime(traj).val(obj.setup.plot_vars(i), :), 'LineStyle', '--', 'LineWidth', 2, 'Color', category20(2, :));
            estimated_values = [estimated_values; obj.init.X_est_runtime(traj).val(obj.setup.plot_vars(i), :)];
        end    

        % Labels
        set(gca, 'fontsize', fontsize);
        ylabel(['x_', num2str(obj.setup.plot_vars(i))]);
        
        xlabel('time[s]');
        legend('True', 'Est');


        % Second subplot: Error
        subplot(2, 1, 2);
        hold on;
        grid on;
        box on;

        for traj = 1:obj.setup.Ntraj
            error = obj.init.X(traj).val(obj.setup.plot_vars(i), :) - obj.init.X_est_runtime(traj).val(obj.setup.plot_vars(i), :);
            plot(obj.setup.time, error, 'LineStyle', '--', 'LineWidth', 2, 'Color', category20(4, :));
        end

        % Calculate Metrics

        if i==1
            errors = true_values(obj.init.Nw_Nts:end) - estimated_values(obj.init.Nw_Nts:end);
            RMSE = sqrt(mean(errors.^2));
            R2 = 1 - sum(errors.^2) / sum((true_values(obj.init.Nw_Nts:end) - mean(true_values(obj.init.Nw_Nts:end))).^2);
            MAPE = mean(abs(errors ./ true_values(obj.init.Nw_Nts:end))) * 100;
    
            % Display Metrics
            metrics_text = {
                ['SOC RMSE: ', num2str(RMSE)];
                ['SOC R²: ', num2str(R2)];
                ['SOC MAPE: ', num2str(MAPE), '%']
            };
            % dim = [.15 .5 .3 .3];
            % annotation('textbox', dim, 'String', metrics_text, 'FitBoxToText', 'on', 'BackgroundColor', 'white');
            disp(metrics_text)
        end

        % Labels
        set(gca, 'fontsize', fontsize);
        ylabel(['Error in x_', num2str(obj.setup.plot_vars(i))]);
        xlabel('time [s]');
        title('Estimation Error');
     

    end    
    
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
            sgtitle('Parameters estimation')
    
            % Plot each parameter in the group
            for i = 1:length(groups{g})
                subplot(length(groups{g}), 1, i);
                hold on
                grid on
                box on
    
                for traj = 1:obj.setup.Ntraj
                    % if 1 || strcmp(obj.setup.DataType, 'simulated')
                    %     plot(obj.setup.time, obj.init.X(traj).val(groups{g}(i), :), 'b--', 'LineWidth', 2);
                    % end
                    plot(obj.setup.time, obj.init.X_est_runtime(traj).val(groups{g}(i), :), 'LineStyle', '--', 'LineWidth', 2, 'Color', category20(2, :));
                end
    
                % labelsval
                set(gca, 'fontsize', fontsize)
                ylabel(ylabels{g}(i,:), 'Interpreter', 'latex');
                end
    
            % if strcmp(obj.setup.DataType, 'simulated')
            %     legend('True', 'Est')
            % else
            %     legend('Stored', 'Est', 'Runtime')
            % end
            xlabel(['time [s]'])
        end
    end

    %%%% plot parameters estimation %%%
    % if ~isempty(obj.setup.plot_params)
    % 
    %     fig_count = fig_count+1;
    %     figure(fig_count)
    %     sgtitle('Parameters estimation')
    %     for i=1:length(obj.setup.plot_params)
    %         subplot(length(obj.setup.plot_params),1,i);
    %         hold on
    %         grid on
    %         box on
    % 
    %         for traj=1:obj.setup.Ntraj
    %             if 1 || strcmp(obj.setup.DataType,'simulated')
    %                 plot(obj.setup.time,obj.init.X(traj).val(obj.setup.plot_params(i),:),'b--','LineWidth',2);
    %             end
    %             plot(obj.setup.time,obj.init.X_est_runtime(traj).val(obj.setup.plot_params(i),:),'g--','LineWidth',2);                                                      
    %         end
    % 
    %         % labels
    %         set(gca,'fontsize', fontsize)             
    %         ylabel(['x_',num2str(obj.setup.plot_params(i))])
    %     end
    % 
    %     if strcat(obj.setup.DataType,'simulated')
    %         legend('True','Est')
    %     else
    %         legend('Stored','Est','Runtime')
    %     end
    %     xlabel(['time [s]'])
    % end
    
%     %%%% plot state estimation error %%%
%     if strcmp(obj.setup.DataType,'simulated')                
%         fig_count = fig_count+1;
%         figure(fig_count)
%         sgtitle('Estimation error - components')
% 
%         for i=1:length(obj.setup.plot_vars)
%             subplot(length(obj.setup.plot_vars),1,i);
%             hold on
%             grid on
%             box on
% 
%             % plot
%             est_error = obj.init.X(1).val(obj.setup.plot_vars(i),:) - obj.init.X_est_runtime(1).val(obj.setup.plot_vars(i),:);
% 
%             log_flag = 1;
%             if ~log_flag
%                 plot(obj.setup.time,est_error,'k','LineWidth',2);
%             else
%                 % log 
% %                     set(gca, 'XScale', 'log')
%                 set(gca, 'YScale', 'log')
%                 plot(obj.setup.time,abs(est_error),'k','LineWidth',2);
%             end            
% 
%             set(gca,'fontsize', fontsize)
%             ylabel(['\delta x_',num2str(obj.setup.plot_vars(i))])
%         end        
%         xlabel('time [s]')        
%     end
% 
%     %%%% plot parameters estimation error %%%
%     if 1 || strcmp(obj.setup.DataType,'simulated')
%         if ~isempty(obj.setup.plot_params)                    
%             fig_count = fig_count+1;
%             figure(fig_count)
%             sgtitle('Estimation error - parameters')
% 
%             for i=1:length(obj.setup.plot_params)
%                 subplot(length(obj.setup.plot_params),1,i);
%                 hold on
%                 grid on
%                 box on
% 
%                 % plot
%                 est_error = obj.init.X(1).val(obj.setup.plot_params(i),:) - obj.init.X_est_runtime(1).val(obj.setup.plot_params(i),:);
% 
%                 log_flag = 1;
%                 if ~log_flag
%                     plot(obj.setup.time,est_error,'b','LineWidth',2);
%                 else
%                     % log 
% %                     set(gca, 'XScale', 'log')
%                     set(gca, 'YScale', 'log')
%                     plot(obj.setup.time,abs(est_error),'b','LineWidth',2);
%                 end
% 
%                 set(gca,'fontsize', fontsize)                
%                 ylabel(['\delta x_',num2str(obj.setup.plot_params(i))])
%             end
%             xlabel('time [s]')
%         end
%     end
    
    %%%% plot state estimation error - norm%%%
    if strcmp(obj.setup.DataType,'simulated')                
        fig_count = fig_count+1;
        figure(fig_count)
        sgtitle('Estimation error state - norm')
        hold on
        grid on
        box on

        % plot
        for iter=1:obj.setup.Niter
            est_error_norm(iter) = norm(obj.init.X(1).val(obj.setup.plot_vars,iter) - obj.init.X_est_runtime(1).val(obj.setup.plot_vars,iter));
        end

        log_flag = 0;
        if ~log_flag
            plot(obj.setup.time,est_error_norm,'k','LineWidth',2);
        else
            % log 
%                     set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')
            plot(obj.setup.time,abs(est_error_norm),'r--','LineWidth',2);
        end

        set(gca,'fontsize', fontsize)
        xlabel('time [s]')
        ylabel('\delta x_norm') 
    end
    
    %%%% plot params estimation error - norm%%%
    if 1 || strcmp(obj.setup.DataType,'simulated')                
        fig_count = fig_count+1;
        figure(fig_count)
        sgtitle('Estimation error params - norm')
        hold on
        grid on
        box on

        % plot
        for iter=1:obj.setup.Niter
            est_error_norm(iter) = norm(obj.init.X(1).val(obj.setup.plot_params,iter) - obj.init.X_est_runtime(1).val(obj.setup.plot_params,iter));
        end

        log_flag = 0;
        if ~log_flag
            plot(obj.setup.time,est_error_norm,'r','LineWidth',2);
        else
            % log 
%                     set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')
            plot(obj.setup.time,abs(est_error_norm),'b--','LineWidth',2);
        end

        set(gca,'fontsize', fontsize)
        xlabel('time [s]')
        ylabel('\delta x_norm') 
    end

    %%%% plot filters %%%%%            
    fig_count = fig_count+1;
    figure(fig_count)
    sgtitle('Filters on measures')            
    ax = zeros(1,3);
    for k=1:obj.setup.J_nterm
        
        % number fo subplots depending on the Nterm
        n_subplot = obj.setup.J_nterm;
        
        % indicize axes
        ax_index = k;
        ax(ax_index)=subplot(n_subplot,1,ax_index);                
        
        % plot
        hold on
        grid on
        
        for traj=1:obj.setup.Ntraj
            for dim=obj.setup.dim_out_compare
                y_plot = obj.setup.J_temp_scale(k)*reshape(obj.init.Y_full_story(traj).val(k,dim,:),size(obj.setup.time));
                if strcmp(obj.setup.DataType,'simulated')
                    ytrue_plot = obj.setup.J_temp_scale(k)*reshape(obj.init.Ytrue_full_story(traj).val(k,dim,:),size(obj.setup.time));
                end
                yhat_plot = obj.setup.J_temp_scale(k)*reshape(obj.init.Yhat_full_story(traj).val(k,dim,:),size(obj.setup.time));

                if strcmp(obj.setup.DataType,'simulated')
                    plot(obj.setup.time,y_plot,'b--');
                end
                plot(obj.setup.time,yhat_plot,'LineStyle', ':', 'LineWidth', 2, 'Color', category20(2, :));
                plot(obj.setup.time,y_plot,'LineStyle', '--', 'LineWidth', 2, 'Color', 'blue');                            

            end                        
        end    
        legend('est','meas')
        set(gca,'fontsize', fontsize)
        ylabel('voltage [V]');
        xlabel('time [s]');

        errors = y_plot(obj.init.Nw_Nts:end) - yhat_plot(obj.init.Nw_Nts:end);
        RMSE = sqrt(mean(errors.^2));
        R2 = 1 - sum(errors.^2) / sum((y_plot(obj.init.Nw_Nts:end) - mean(y_plot(obj.init.Nw_Nts:end))).^2);
        MAPE = mean(abs(errors ./ y_plot(obj.init.Nw_Nts:end))) * 100;

        % Display Metrics
        metrics_text = {
            ['Voltage RMSE: ', num2str(RMSE)];
            ['Voltage R²: ', num2str(R2)];
            ['Voltage MAPE: ', num2str(MAPE), '%']
        };
        % dim = [.15 .5 .3 .3];
        % annotation('textbox', dim, 'String', metrics_text, 'FitBoxToText', 'on', 'BackgroundColor', 'white');
        disp(metrics_text)
        
    end
    linkaxes(ax,'x');  
    
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
    fig_count = fig_count+1;
    figure(fig_count)
    cwt(y_meas,obj.init.wvname,1/obj.setup.Ts,'VoicesPerOctave',obj.init.Nv,'FrequencyLimits',obj.init.FLIMITS);
            
    
end