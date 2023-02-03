function [Valor_compensacao] = calc_compensacao_consumidor2(van, vbn, vcn)
    
    % Parametrizacao
    DRP_madrugada = 0;
    DRP_manha = 0;
    DRP_tarde = 0;
    DRP_noite = 0;
    DRPm = 0;
    DRC_madrugada = 0;
    DRC_manha = 0;
    DRC_tarde = 0;
    DRC_noite = 0;
    DRCm = 0;
    EUSD = 10;
   
    % Seleciona valores uteis   
    madrugada = [van(1), vbn(1), vcn(1)];
    min_madrugada = min(madrugada);
    max_madrugada = max(madrugada);

    manha = [van(2), vbn(2), vcn(2)];
    min_manha = min(manha);
    max_manha = max(manha);

    tarde = [van(3), vbn(3), vcn(3)];
    min_tarde = min(tarde);
    max_tarde = max(tarde);

    noite = [van(4), vbn(4), vcn(4)];
    min_noite = min(noite);
    max_noite = max(noite);


    % ###Madrugada###
    % DRC
    if min_madrugada < 0.9
        DRC_madrugada = DRC_madrugada + 25;
    elseif max_madrugada > 1.05
        DRC_madrugada = DRC_madrugada + 25;
    % DRP
    elseif min_madrugada < 0.93
        DRP_madrugada = DRP_madrugada + 25;
    end

    % k1
    if DRP_madrugada <= DRPm
        k1_madrugada= 0;
    else
        k1_madrugada = 3;
    end
        
    % k2
    if DRC_madrugada <= DRCm
        k2_madrugada = 0;
    else
        k2_madrugada = 5;
    end

    % ###Manha###
    % DRC
    if min_manha < 0.9
        DRC_manha = DRC_manha + 25;
    elseif max_manha > 1.05
        DRC_manha = DRC_manha + 25;
    % DRP
    elseif min_manha < 0.93
        DRP_manha = DRP_manha + 25;
    end

    % k1
    if DRP_manha <= DRPm
        k1_manha = 0;
    else
        k1_manha = 3;
    end
        
    % k2
    if DRC_manha <= DRCm
        k2_manha = 0;
    else
        k2_manha = 5;
    end

    % ###Tarde###
    % DRC
    if min_tarde < 0.9
        DRC_tarde = DRC_tarde + 25;
    elseif max_tarde > 1.05
        DRC_tarde = DRC_tarde + 25;
    % DRP
    elseif min_tarde < 0.93
        DRP_tarde = DRP_tarde + 25;
    end

    % k1
    if DRP_tarde <= DRPm
        k1_tarde = 0;
    else
        k1_tarde = 3;
    end
        
    % k2
    if DRC_tarde <= DRCm
        k2_tarde = 0;
    else
        k2_tarde = 5;
    end

    % ###Noite###
    % DRC
    if min_noite < 0.9
        DRC_noite = DRC_noite + 25;
    elseif max_noite > 1.05
        DRC_noite = DRC_noite + 25;
    % DRP
    elseif min_noite < 0.93
        DRP_noite = DRP_noite + 25;
    end

    % k1
    if DRP_noite <= DRPm
        k1_noite = 0;
    else
        k1_noite = 3;
    end
        
    % k2
    if DRC_noite <= DRCm
        k2_noite = 0;
    else
        k2_noite = 5;
    end

    % Valor
    Valor_madrugada = (k1_madrugada * ((DRP_madrugada - DRPm) / 100) + k2_madrugada * ((DRC_madrugada - DRP_madrugada) / 100)) * EUSD;
    Valor_manha = (k1_manha * ((DRP_manha - DRPm) / 100) + k2_manha * ((DRC_manha - DRP_manha) / 100)) * EUSD;
    Valor_tarde = (k1_tarde * ((DRP_tarde - DRPm) / 100) + k2_tarde * ((DRC_tarde - DRP_tarde) / 100)) * EUSD;
    Valor_noite = (k1_noite * ((DRP_noite - DRPm) / 100) + k2_noite * ((DRC_noite - DRP_noite) / 100)) * EUSD;
    
    Valor_compensacao = [Valor_madrugada, Valor_manha, Valor_tarde, Valor_noite];

end
