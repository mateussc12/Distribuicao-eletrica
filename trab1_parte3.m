clc; clear;
% Parametrizacao
z1 = 0.2053 + 1j*0.3753; % Ohm/km
z2 = 2 * z1; % Ohm/km
tolerancia = 0.000001;
V_se = 1; % pu

% Figura 1
d12 = 0.5; % km
d23 = 0.5; % km
d24 = 2; % km
d35 = 0.75; % km
d36 = 1.5; % km

% Tabela 4
S2_EA_ER = [800e6 + 1j*600e6, 1200e6 + 1j*400e6, 1500e6 + 1j*1000e6] / (24 * 30);   % MWh/h
S4_EA_ER = [1500e6 + 1j*600e6, 1800e6 + 1j*1200e6, 2000e6 + 1j*1200e6] / (24 * 30); % MWh/h
S5_EA_ER = [2500e6 + 1j*1000e6, 3000e6 + 1j*1500e6, 1000e6 + 1j*800e6] / (24 * 30); % MWh/h
S6_EA_ER = [2500e6 + 1j*500e6, 1500e6 + 1j*1500e6, 1800e6 + 1j*1400e6] / (24 * 30); % MWh/h

% Tabela 5
DA = [0.6, 0.9, 1, 0.8];
DR = [0.4, 0.5, 0.6, 0.6];

% Tabela 6
zip_S2 = [0.2, 0.3, 0.5];
zip_S4 = [0.2, 0.6, 0.2];
zip_S5 = [1, 0, 0];
zip_S6 = [0, 0.5, 0.5];

% Valores de base
S_base_tri = 10e6; % VA
S_base_mono = S_base_tri / 3; % VA
V_base_linha = 11.4e3; % V
V_base_fase = V_base_linha / sqrt(3); % V
z1_pu = (S_base_mono / V_base_fase^2) * z1; % pu/km
z2_pu = (S_base_mono / V_base_fase^2) * z2; % pu/km

% Valores base
for i = 1:length(DA)
S2_S_fases_abc(i, :) = real(S2_EA_ER) * DA(i) / S_base_mono + 1j * imag(S2_EA_ER) * DR(i) / S_base_mono; % pu
S4_S_fases_abc(i, :) = real(S4_EA_ER) * DA(i) / S_base_mono + 1j * imag(S4_EA_ER) * DR(i) / S_base_mono; % pu
S5_S_fases_abc(i, :) = real(S5_EA_ER) * DA(i) / S_base_mono + 1j * imag(S5_EA_ER) * DR(i) / S_base_mono; % pu
S6_S_fases_abc(i, :) = real(S6_EA_ER) * DA(i) / S_base_mono + 1j * imag(S6_EA_ER) * DR(i) / S_base_mono; % pu
end

% Questao 1
v_a_e_pu = [1, 1, 1, 1]; % pu
v_b_e_pu = [-0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2]; % pu
v_c_e_pu = [-0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2]; % pu
v_n_e_pu = [0, 0, 0, 0]; % pu

v_a_s_pu = [1, 1, 1, 1]; % pu
v_b_s_pu = [-0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2]; % pu
v_c_s_pu = [-0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2]; % pu
v_n_s_pu = [0, 0, 0, 0]; % pu

% Inicializando as tensões nas barras
v1a = v_a_s_pu; v1b = v_b_s_pu; v1c = v_c_s_pu; v1n = v_n_s_pu;
v2a = v_a_s_pu; v2b = v_b_s_pu; v2c = v_c_s_pu; v2n = v_n_s_pu;
v3a = v_a_s_pu; v3b = v_b_s_pu; v3c = v_c_s_pu; v3n = v_n_s_pu;
v4a = v_a_s_pu; v4b = v_b_s_pu; v4c = v_c_s_pu; v4n = v_n_s_pu;
v5a = v_a_s_pu; v5b = v_b_s_pu; v5c = v_c_s_pu; v5n = v_n_s_pu;
v6a = v_a_s_pu; v6b = v_b_s_pu; v6c = v_c_s_pu; v6n = v_n_s_pu;

v_ref = 1; % pu

erro_2a = [1, 1, 1, 1];
erro_2b = [1, 1, 1, 1];
erro_2c = [1, 1, 1, 1];
erro_2n = [1, 1, 1, 1];

erro_4a = [1, 1, 1, 1];
erro_4b = [1, 1, 1, 1];
erro_4c = [1, 1, 1, 1];
erro_4n = [1, 1, 1, 1];

erro_5a = [1, 1, 1, 1];
erro_5b = [1, 1, 1, 1];
erro_5c = [1, 1, 1, 1];
erro_5n = [1, 1, 1, 1];

erro_6a = [1, 1, 1, 1];
erro_6b = [1, 1, 1, 1];
erro_6c = [1, 1, 1, 1];
erro_6n = [1, 1, 1, 1];

iter = 1;
% Calcula correntes nas barras   
while (((erro_2a(iter, :) > tolerancia) | (erro_2b(iter, :) > tolerancia) | (erro_2c(iter, :) > tolerancia) | (erro_2n(iter, :) > tolerancia)) | ((erro_4a(iter, :) > tolerancia) | (erro_4b(iter, :) > tolerancia) | (erro_4c(iter, :) > tolerancia) | (erro_4n(iter, :) > tolerancia)) | ((erro_5a(iter, :) > tolerancia) | (erro_5b(iter, :) > tolerancia) | (erro_5c(iter, :) > tolerancia) | (erro_5n(iter, :) > tolerancia)) | ((erro_6a(iter, :) > tolerancia) | (erro_6b(iter, :) > tolerancia) | (erro_6c(iter, :) > tolerancia) | (erro_6n(iter, :) > tolerancia)))
    
    [S2_fase_a, S2_fase_b, S2_fase_c] = ordem_potencia(S2_S_fases_abc);
    i2a(iter, :) = (zip_S2(3) .* (conj(S2_fase_a) ./ conj(v2a(iter, :) - v2n(iter, :)))) + (zip_S2(2) .* (conj(S2_fase_a) ./ conj(v_ref)) .* ((v2a(iter, :) - v2n(iter, :)) ./ abs((v2a(iter, :) - v2n(iter, :))))) + (zip_S2(1) .* (conj(S2_fase_a)./((v_ref)^2)) .* (v2a(iter, :) - v2n(iter, :)));
    i2b(iter, :) = (zip_S2(3) .* (conj(S2_fase_b) ./ conj(v2b(iter, :) - v2n(iter, :)))) + (zip_S2(2) .* (conj(S2_fase_b) ./ conj(v_ref)) .* ((v2b(iter, :) - v2n(iter, :)) ./ abs((v2b(iter, :) - v2n(iter, :))))) + (zip_S2(1) .* (conj(S2_fase_b)./((v_ref)^2)) .* (v2b(iter, :) - v2n(iter, :)));
    i2c(iter, :) = (zip_S2(3) .* (conj(S2_fase_c) ./ conj(v2c(iter, :) - v2n(iter, :)))) + (zip_S2(2) .* (conj(S2_fase_c) ./ conj(v_ref)) .* ((v2c(iter, :) - v2n(iter, :)) ./ abs((v2c(iter, :) - v2n(iter, :))))) + (zip_S2(1) .* (conj(S2_fase_c)./((v_ref)^2)) .* (v2c(iter, :) - v2n(iter, :)));
    i2n(iter, :) = - (i2a(iter, :) + i2b(iter, :) + i2c(iter, :));

    [S4_fase_a, S4_fase_b, S4_fase_c] = ordem_potencia(S4_S_fases_abc);
    i4a(iter, :) = (zip_S4(3) .* (conj(S4_fase_a) ./ conj(v4a(iter, :) - v4n(iter, :)))) + (zip_S4(2) .* (conj(S4_fase_a) ./ conj(v_ref)) .* ((v4a(iter, :) - v4n(iter, :)) ./ abs((v4a(iter, :) - v4n(iter, :))))) + (zip_S4(1) .* (conj(S4_fase_a)./((v_ref)^2)) .* (v4a(iter, :) - v4n(iter, :)));
    i4b(iter, :) = (zip_S4(3) .* (conj(S4_fase_b) ./ conj(v4b(iter, :) - v4n(iter, :)))) + (zip_S4(2) .* (conj(S4_fase_b) ./ conj(v_ref)) .* ((v4b(iter, :) - v4n(iter, :)) ./ abs((v4b(iter, :) - v4n(iter, :))))) + (zip_S4(1) .* (conj(S4_fase_b)./((v_ref)^2)) .* (v4b(iter, :) - v4n(iter, :)));
    i4c(iter, :) = (zip_S4(3) .* (conj(S4_fase_c) ./ conj(v4c(iter, :) - v4n(iter, :)))) + (zip_S4(2) .* (conj(S4_fase_c) ./ conj(v_ref)) .* ((v4c(iter, :) - v4n(iter, :)) ./ abs((v4c(iter, :) - v4n(iter, :))))) + (zip_S4(1) .* (conj(S4_fase_c)./((v_ref)^2)) .* (v4c(iter, :) - v4n(iter, :)));
    i4n(iter, :) = - (i4a(iter, :) + i4b(iter, :) + i4c(iter, :));

    [S5_fase_a, S5_fase_b, S5_fase_c] = ordem_potencia(S5_S_fases_abc);
    i5a(iter, :) = (zip_S5(3) .* (conj(S5_fase_a) ./ conj(v5a(iter, :) - v5n(iter, :)))) + (zip_S5(2) .* (conj(S5_fase_a) ./ conj(v_ref)) .* ((v5a(iter, :) - v5n(iter, :)) ./ abs((v5a(iter, :) - v5n(iter, :))))) + (zip_S5(1) .* (conj(S5_fase_a)./((v_ref)^2)) .* (v5a(iter, :) - v5n(iter, :)));
    i5b(iter, :) = (zip_S5(3) .* (conj(S5_fase_b) ./ conj(v5b(iter, :) - v5n(iter, :)))) + (zip_S5(2) .* (conj(S5_fase_b) ./ conj(v_ref)) .* ((v5b(iter, :) - v5n(iter, :)) ./ abs((v5b(iter, :) - v5n(iter, :))))) + (zip_S5(1) .* (conj(S5_fase_b)./((v_ref)^2)) .* (v5b(iter, :) - v5n(iter, :)));
    i5c(iter, :) = (zip_S5(3) .* (conj(S5_fase_c) ./ conj(v5c(iter, :) - v5n(iter, :)))) + (zip_S5(2) .* (conj(S5_fase_c) ./ conj(v_ref)) .* ((v5c(iter, :) - v5n(iter, :)) ./ abs((v5c(iter, :) - v5n(iter, :))))) + (zip_S5(1) .* (conj(S5_fase_c)./((v_ref)^2)) .* (v5c(iter, :) - v5n(iter, :)));
    i5n(iter, :) = - (i5a(iter, :) + i5b(iter, :) + i5c(iter, :));

    [S6_fase_a, S6_fase_b, S6_fase_c] = ordem_potencia(S6_S_fases_abc);
    i6a(iter, :) = (zip_S6(3) .* (conj(S6_fase_a) ./ conj(v6a(iter, :) - v6n(iter, :)))) + (zip_S6(2) .* (conj(S6_fase_a) ./ conj(v_ref)) .* ((v6a(iter, :) - v6n(iter, :)) ./ abs((v6a(iter, :) - v6n(iter, :))))) + (zip_S6(1) .* (conj(S6_fase_a)./((v_ref)^2)) .* (v6a(iter, :) - v6n(iter, :)));
    i6b(iter, :) = (zip_S6(3) .* (conj(S6_fase_b) ./ conj(v6b(iter, :) - v6n(iter, :)))) + (zip_S6(2) .* (conj(S6_fase_b) ./ conj(v_ref)) .* ((v6b(iter, :) - v6n(iter, :)) ./ abs((v6b(iter, :) - v6n(iter, :))))) + (zip_S6(1) .* (conj(S6_fase_b)./((v_ref)^2)) .* (v6b(iter, :) - v6n(iter, :)));
    i6c(iter, :) = (zip_S6(3) .* (conj(S6_fase_c) ./ conj(v6c(iter, :) - v6n(iter, :)))) + (zip_S6(2) .* (conj(S6_fase_c) ./ conj(v_ref)) .* ((v6c(iter, :) - v6n(iter, :)) ./ abs((v6c(iter, :) - v6n(iter, :))))) + (zip_S6(1) .* (conj(S6_fase_c)./((v_ref)^2)) .* (v6c(iter, :) - v6n(iter, :)));
    i6n(iter, :) = - (i6a(iter, :) + i6b(iter, :) + i6c(iter, :));
    
    % Calcula correntes nos trechos
    i2_3a(iter, :) = i5a(iter, :) + i6a(iter, :);
    i2_3b(iter, :) = i5b(iter, :) + i6b(iter, :);
    i2_3c(iter, :) = i5c(iter, :) + i6c(iter, :);
    i2_3n(iter, :) = i5n(iter, :) + i6n(iter, :);

    i1_2a(iter, :) = i2_3a(iter, :) +  i2a(iter, :) + i4a(iter, :);
    i1_2b(iter, :) = i2_3b(iter, :) +  i2b(iter, :) + i4b(iter, :);
    i1_2c(iter, :) = i2_3c(iter, :) +  i2c(iter, :) + i4c(iter, :);
    i1_2n(iter, :) = i2_3n(iter, :) +  i2n(iter, :) + i4n(iter, :);

    i2_4a(iter, :) = i4a(iter, :);
    i2_4b(iter, :) = i4b(iter, :);
    i2_4c(iter, :) = i4c(iter, :);
    i2_4n(iter, :) = i4n(iter, :);

    i3_5a(iter, :) = i5a(iter, :);
    i3_5b(iter, :) = i5b(iter, :);
    i3_5c(iter, :) = i5c(iter, :);
    i3_5n(iter, :) = i5n(iter, :);

    i3_6a(iter, :) = i6a(iter, :);
    i3_6b(iter, :) = i6b(iter, :);
    i3_6c(iter, :) = i6c(iter, :);
    i3_6n(iter, :) = i6n(iter, :);

    % Calcula a queda de tensão nos trechos
    delta_1_2a(iter, :) = d12 * z1_pu * i1_2a(iter, :);
    delta_1_2b(iter, :) = d12 * z1_pu * i1_2b(iter, :);
    delta_1_2c(iter, :) = d12 * z1_pu * i1_2c(iter, :);
    delta_1_2n(iter, :) = d12 * z1_pu * i1_2n(iter, :);

    delta_2_3a(iter, :) = d23 * z1_pu * i2_3a(iter, :);
    delta_2_3b(iter, :) = d23 * z1_pu * i2_3b(iter, :);
    delta_2_3c(iter, :) = d23 * z1_pu * i2_3c(iter, :);
    delta_2_3n(iter, :) = d23 * z1_pu * i2_3n(iter, :);

    delta_2_4a(iter, :) = d24 * z2_pu * i2_4a(iter, :);
    delta_2_4b(iter, :) = d24 * z2_pu * i2_4b(iter, :);
    delta_2_4c(iter, :) = d24 * z2_pu * i2_4c(iter, :);
    delta_2_4n(iter, :) = d24 * z2_pu * i2_4n(iter, :);

    delta_3_5a(iter, :) = d35 * z2_pu * i3_5a(iter, :);
    delta_3_5b(iter, :) = d35 * z2_pu * i3_5b(iter, :);
    delta_3_5c(iter, :) = d35 * z2_pu * i3_5c(iter, :);
    delta_3_5n(iter, :) = d35 * z2_pu * i3_5n(iter, :);

    delta_3_6a(iter, :) = d36 * z2_pu * i3_6a(iter, :) / 2;
    delta_3_6b(iter, :) = d36 * z2_pu * i3_6b(iter, :) / 2;
    delta_3_6c(iter, :) = d36 * z2_pu * i3_6c(iter, :) / 2;
    delta_3_6n(iter, :) = d36 * z2_pu * i3_6n(iter, :) / 2;

    % Calcula as tensões dos trechos
    v2a(iter + 1, :) = v_a_e_pu - delta_1_2a(iter, :);
    v2b(iter + 1, :) = v_b_e_pu - delta_1_2b(iter, :);
    v2c(iter + 1, :) = v_c_e_pu - delta_1_2c(iter, :);
    v2n(iter + 1, :) = v_n_e_pu - delta_1_2n(iter, :);

    v4a(iter + 1, :) = v_a_e_pu - delta_2_4a(iter, :);
    v4b(iter + 1, :) = v_b_e_pu - delta_2_4b(iter, :);
    v4c(iter + 1, :) = v_c_e_pu - delta_2_4c(iter, :);
    v4n(iter + 1, :) = v_n_e_pu - delta_2_4n(iter, :);

    v5a(iter + 1, :) = v_a_e_pu - delta_3_5a(iter, :);
    v5b(iter + 1, :) = v_b_e_pu - delta_3_5b(iter, :);
    v5c(iter + 1, :) = v_c_e_pu - delta_3_5c(iter, :);
    v5n(iter + 1, :) = v_n_e_pu - delta_3_5n(iter, :);

    v6a(iter + 1, :) = v_a_e_pu - delta_3_6a(iter, :);
    v6b(iter + 1, :) = v_b_e_pu - delta_3_6b(iter, :);
    v6c(iter + 1, :) = v_c_e_pu - delta_3_6c(iter, :);
    v6n(iter + 1, :) = v_n_e_pu - delta_3_6n(iter, :);

    % Calcula os Erros
    erro_2a(iter + 1, :) = abs(v2a(iter + 1, :) - v2a(iter, :));
    erro_2b(iter + 1, :) = abs(v2b(iter + 1, :) - v2b(iter, :));
    erro_2c(iter + 1, :) = abs(v2c(iter + 1, :) - v2c(iter, :));
    erro_2n(iter + 1, :) = abs(v2n(iter + 1, :) - v2n(iter, :));

    erro_4a(iter + 1, :) = abs(v4a(iter + 1, :) - v4a(iter, :));
    erro_4b(iter + 1, :) = abs(v4b(iter + 1, :) - v4b(iter, :));
    erro_4c(iter + 1, :) = abs(v4c(iter + 1, :) - v4c(iter, :));
    erro_4n(iter + 1, :) = abs(v4n(iter + 1, :) - v4n(iter, :));

    erro_5a(iter + 1, :) = abs(v5a(iter + 1, :) - v5a(iter, :));
    erro_5b(iter + 1, :) = abs(v5b(iter + 1, :) - v5b(iter, :));
    erro_5c(iter + 1, :) = abs(v5c(iter + 1, :) - v5c(iter, :));
    erro_5n(iter + 1, :) = abs(v5n(iter + 1, :) - v5n(iter, :));

    erro_6a(iter + 1, :) = abs(v6a(iter + 1, :) - v6a(iter, :));
    erro_6b(iter + 1, :) = abs(v6b(iter + 1, :) - v6b(iter, :));
    erro_6c(iter + 1, :) = abs(v6c(iter + 1, :) - v6c(iter, :));
    erro_6n(iter + 1, :) = abs(v6n(iter + 1, :) - v6n(iter, :));

    iter = iter + 1;

end

% Calcula as tens�es de fase-neutro
v2an = v2a(end, :) - v2n(end, :);
v2bn = v2b(end, :) - v2n(end, :);
v2cn = v2c(end, :) - v2n(end, :);

v4an = v4a(end, :) - v4n(end, :);
v4bn = v4b(end, :) - v4n(end, :);
v4cn = v4c(end, :) - v4n(end, :);

v5an = v5a(end, :) - v5n(end, :);
v5bn = v5b(end, :) - v5n(end, :);
v5cn = v5c(end, :) - v5n(end, :);

v6an = v6a(end, :) - v6n(end, :);
v6bn = v6b(end, :) - v6n(end, :);
v6cn = v6c(end, :) - v6n(end, :);

v2an = [round(abs(v2an), 4), round(rad2deg(angle(v2an)), 2)]
v2bn = [round(abs(v2bn), 4), round(rad2deg(angle(v2cn)), 2)]
v2cn = [round(abs(v2cn), 4), round(rad2deg(angle(v2cn)), 2)]

v4an = [round(abs(v4an), 4), round(rad2deg(angle(v4an)), 2)]
v4bn = [round(abs(v4bn), 4), round(rad2deg(angle(v4bn)), 2)]
v4cn = [round(abs(v4cn), 4), round(rad2deg(angle(v4cn)), 2)]

v5an = [round(abs(v5an), 4), round(rad2deg(angle(v5an)), 2)]
v5bn = [round(abs(v5bn), 4), round(rad2deg(angle(v5bn)), 2)]
v5cn = [round(abs(v5cn), 4), round(rad2deg(angle(v5cn)), 2)]

v6an = [round(abs(v6an), 4), round(rad2deg(angle(v6an)), 2)]
v6bn = [round(abs(v6bn), 4), round(rad2deg(angle(v6bn)), 2)]
v6cn = [round(abs(v6cn), 4), round(rad2deg(angle(v6cn)), 2)]

% Calcula as perdas
pot_entrada = real((v_a_e_pu .* conj(i1_2a(end, :))) + (v_b_e_pu .* conj(i1_2b(end, :))) + (v_c_e_pu .* conj(i1_2c(end, :))))

%pot_saida_col1 = real(sum(S2_S_fases_abc(1, :)) + sum(S4_S_fases_abc(1, :)) + sum(S5_S_fases_abc(1, :)) + sum(S6_S_fases_abc(1, :)));
%pot_saida_col2 = real(sum(S2_S_fases_abc(2, :)) + sum(S4_S_fases_abc(2, :)) + sum(S5_S_fases_abc(2, :)) + sum(S6_S_fases_abc(2, :)));
%pot_saida_col3 = real(sum(S2_S_fases_abc(3, :)) + sum(S4_S_fases_abc(3, :)) + sum(S5_S_fases_abc(3, :)) + sum(S6_S_fases_abc(3, :)));
%pot_saida_col4 = real(sum(S2_S_fases_abc(4, :)) + sum(S4_S_fases_abc(4, :)) + sum(S5_S_fases_abc(4, :)) + sum(S6_S_fases_abc(4, :)));
%pot_saida = [pot_saida_col1, pot_saida_col2, pot_saida_col3, pot_saida_col4]

pot_saida = ( (v2a + v2b + v2c + v2n) * conj(i)

perdas = ((pot_entrada  - pot_saida) ./ pot_saida) * 100
