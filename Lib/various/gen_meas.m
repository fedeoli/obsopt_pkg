%%
function x_out = gen_meas(params,x,Npoint,sigma)    
             
    NSeg = length(params.input_data.SOC)-1;
    NPointsSeg = floor(Npoint/NSeg);
    
    % initialize state
    x_out = zeros(length(x),NSeg*NPointsSeg);
    
    % create cloud
    for i=1:NSeg
        for j=1:NPointsSeg
            % generate SOC        
            x_out(1,(i-1)*NPointsSeg+j) =  unifrnd(params.input_data.SOC(i),params.input_data.SOC(i+1));            
            
            % generate points (noiseless)
            ref(3) = interp1(params.input_data.SOC(i:i+1), params.input_data.OCV(i:i+1), x_out(1,(i-1)*NPointsSeg+j));
            ref(4) = interp1(params.input_data.SOC(i:i+1), params.input_data.R0(i:i+1), x_out(1,(i-1)*NPointsSeg+j));
            ref(5) = interp1(params.input_data.SOC(i:i+1), params.input_data.R1(i:i+1), x_out(1,(i-1)*NPointsSeg+j));
            ref(6) = interp1(params.input_data.SOC(i:i+1), params.input_data.C1(i:i+1), x_out(1,(i-1)*NPointsSeg+j));
            
            x_out(3,(i-1)*NPointsSeg+j) = ref(3)*(1+sigma*randn);
            x_out(4,(i-1)*NPointsSeg+j) = ref(4)*(1+sigma*randn);
            x_out(5,(i-1)*NPointsSeg+j) = ref(5)*(1+sigma*randn);
            x_out(6,(i-1)*NPointsSeg+j) = ref(6)*(1+sigma*randn);
        
        end
    end  

    [x_out(1,:), I] = sort(x_out(1,:));
    x_out(3,:) = x_out(3,I);
    x_out(4,:) = x_out(4,I);
    x_out(5,:) = x_out(5,I);
    x_out(6,:) = x_out(6,I);
end