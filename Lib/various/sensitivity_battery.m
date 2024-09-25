%%

function out = sensitivity_battery(params)

    %% generate data
    sigma_points = 0;
    Npoint = 100;
    x_start = zeros(6,1);
    x_out = gen_meas(params,x_start,Npoint,sigma_points);
    SOC = x_out(1,:);
    Voc = x_out(3,:);
    R0 = x_out(4,:);
    R1 = x_out(5,:);
    C1 = x_out(6,:);
    
    %% which window I use?
    % step = 60;
    % Istart = randi(length(SOC)-step,1);    
    % Iend = Istart+step;
    Istart = 1;
    Iend = length(SOC);
    window = Istart:Iend;
    
    %% get initial estimate
    Npoly = 4;    
    
    % coefficients
    % Voc_P = flip(obs.init.X_est_runtime.val(update_vars(1,:),1));
    % R0_P =  flip(obs.init.X_est_runtime.val(update_vars(2,:),1));
    % R1_P =  flip(obs.init.X_est_runtime.val(update_vars(3,:),1));
    % C1_P =  flip(obs.init.X_est_runtime.val(update_vars(4,:),1));
    Voc_P = polyfit(SOC(window),Voc(window),Npoly);
    R0_P =  polyfit(SOC(window),R0(window),Npoly);
    R1_P =  polyfit(SOC(window),R1(window),Npoly);
    C1_P =  polyfit(SOC(window),C1(window),Npoly);
    
    
    % estimation
    Voc_est = polyval(Voc_P,SOC(window));
    R0_est =  polyval(R0_P,SOC(window));
    R1_est =  polyval(R1_P,SOC(window));
    C1_est =  polyval(C1_P,SOC(window));
    
    %% check effect of parameters
    Ncheck = 10;
    Sigma = 1*1e-2;
    % param_check_id = [1:Npoly];
    param_check_id = 1;
    Nparam_check = length(param_check_id);
    
    for i=1:Ncheck
    
        % init params
        Voc_P_check(i,:) = Voc_P;
        R0_P_check(i,:) = R0_P;
        R1_P_check(i,:) = R1_P;
        C1_P_check(i,:) = C1_P;
    
        % create disturb
        Voc_P_check(i,param_check_id) = Voc_P(param_check_id).*(1 + Sigma*randn([1 Nparam_check]));
        R0_P_check(i,param_check_id) = R0_P(param_check_id).*(1 + Sigma*randn([1 Nparam_check]));
        R1_P_check(i,param_check_id) = R1_P(param_check_id).*(1 + Sigma*randn([1 Nparam_check]));
        C1_P_check(i,param_check_id) = C1_P(param_check_id).*(1 + Sigma*randn([1 Nparam_check]));
    
        % estimate - normal
        Voc_est_check (i,:) = polyval(Voc_P_check(i,:),SOC(window));
        R0_est_check(i,:)   =  polyval(R0_P_check(i,:),SOC(window));
        R1_est_check(i,:)   =  polyval(R1_P_check(i,:),SOC(window));
        C1_est_check(i,:)   =  polyval(C1_P_check(i,:),SOC(window)); 

        % error
        out.E_Voc(:,i) = (Voc_est - Voc_est_check(i,:)).^2;
        out.E_R0(:,i) = (R0_est - R0_est_check(i,:)).^2;
        out.E_R1(:,i) = (R1_est - R1_est_check(i,:)).^2;
        out.E_C1(:,i) = (C1_est - C1_est_check(i,:)).^2;

        % RMSE
        out.E_Voc_RMSE(i) = sqrt(mean(out.E_Voc(:,i)));
        out.E_R0_RMSE(i) = sqrt(mean(out.E_R0(:,i)));
        out.E_R1_RMSE(i) = sqrt(mean(out.E_R1(:,i)));
        out.E_C1_RMSE(i) = sqrt(mean(out.E_C1(:,i)));
    
    end
    
    %% plots
    
    % Voc
    figure(1)
    title('Voc')
    xlabel('SOC')
    ylabel('Voc')
    hold on
    
    plot(SOC(window),Voc(window),'bo')
    plot(SOC(window),Voc_est,'r--')
    for i=1:Ncheck
        plot(SOC(window),Voc_est_check(i,:),'k:')
    end
    
    % R0
    figure(2)
    title('R0')
    xlabel('SOC')
    ylabel('R0')
    hold on
    
    plot(SOC(window),R0(window),'bo')
    plot(SOC(window),R0_est,'r--')
    for i=1:Ncheck
        plot(SOC(window),R0_est_check(i,:),'k:')
    end
    
    % R1
    figure(3)
    title('R1')
    xlabel('SOC')
    ylabel('R1')
    hold on
    
    plot(SOC(window),R1(window),'bo')
    plot(SOC(window),R1_est,'r--')
    for i=1:Ncheck
        plot(SOC(window),R1_est_check(i,:),'k:')
    end
    
    % R1
    figure(4)
    title('C1')
    xlabel('SOC')
    ylabel('C1')
    hold on
    
    plot(SOC(window),C1(window),'bo')
    plot(SOC(window),C1_est,'r--')
    for i=1:Ncheck
        plot(SOC(window),C1_est_check(i,:),'k:')
    end
end