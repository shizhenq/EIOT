function S_E=EIOT_MODEL_UPDATE(X_CAL,Y_CAL,X_AUG,Y_AUG,S_hat)

    
    [U,S,V]=svd(X_CAL-Y_CAL*S_hat);
    S_thres=S(1,1);
    
    X_CAL_NEW=[X_CAL; X_AUG];
    Y_CAL_NEW=[Y_CAL;Y_AUG];

    [U_NEW,S_NEW,V_NEW]=svd(X_CAL_NEW-Y_CAL_NEW*S_hat);
    SR=X_CAL_NEW-Y_CAL_NEW*S_hat;
    temp_S=diag(S_NEW(1:size(X_CAL_NEW,1),1:size(X_CAL_NEW,1)));    
%     [a,b]=find(abs(diff(temp_S))>=S_thres);
    [a,b]=find(temp_S>=S_thres);
    if isempty (a)
        S_E=S_hat;
        SR_NEW=SR;
    else
        
        S_I = V_NEW(:,1:length(a));
        S_short = S_NEW(1:length(a),1:length(a));
        r_I     = U_NEW(:,1:length(a))*S_short;
        S_E = [S_hat;S_I'];
        
        SR_NEW  = X_CAL_NEW-[Y_CAL_NEW r_I]*S_E;
    end

    H=figure;
    figure(H)
    subplot(1,2,1)
    plot(1:20,diag(S(1:20,1:20)),'b.-')
    hold on
    plot(1:20,temp_S(1:20),'ro-')
    
    subplot(1,2,2)
    plot(1:size(X_CAL,2),SR,'b-')
    hold on
    plot(1:size(X_CAL,2),SR_NEW,'r.')
    



