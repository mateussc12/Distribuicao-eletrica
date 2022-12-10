clc; clear;
% Parametrizacao
Vse = 11.4e3 / sqrt(3); % V
s_base = 100; % VA
z = 0.2053 + 1j*0.3753; % Ohm/km
d = 5; % km
tolerancia = 0.000001;
erro_a = [1, 1, 1, 1];
erro_b = [1, 1, 1, 1];
erro_c = [1, 1, 1, 1];
erro_n = [1, 1, 1, 1];

% 1° Tabela
EA_mes_fase_a = 1000e6; % Wh/mes
EA_mes_fase_b = 400e6; % Wh/mes
EA_mes_fase_c = 1200e6; % Wh/mes
ER_mes_fase_a = 700e6; % Wh/mes
ER_mes_fase_b = 400e6; % Wh/mes
ER_mes_fase_c = 900e6; % Wh/mes

% 2° Tabela
DA = [0.5, 0.8, 1.3, 0.7];
DR = [0.4, 0.6, 1, 0.5];

% 3° tabela
kz = 0.5;
ki = 0.2;
kp = 0.3;

% Valores de base
S_base_tri = 10e6; % VA
S_base_mono = S_base_tri / 3; % VA
V_base_linha = 11.4e3; % V
V_base_fase = V_base_linha / sqrt(3); % V

% Demanda ativa e reativa do consumidor para os 4 periodos
EA_dia_fase_a = EA_mes_fase_a / (24 * 30); % Wh/h
ER_dia_fase_a = ER_mes_fase_a / (24 * 30); % Wh/h
EA_dia_fase_b = EA_mes_fase_b / (24 * 30); % Wh/h
ER_dia_fase_b = ER_mes_fase_b / (24 * 30); % Wh/h
EA_dia_fase_c = EA_mes_fase_c / (24 * 30); % Wh/h
ER_dia_fase_c = ER_mes_fase_c / (24 * 30); % Wh/h

% Calculos das demanada por fase dos intervalos de tempo 
DA_fase_a = EA_dia_fase_a * DA / S_base_mono; % pu
DR_fase_a = ER_dia_fase_a * DR / S_base_mono; % pu

DA_fase_b = EA_dia_fase_b * DA / S_base_mono; % pu
DR_fase_b = ER_dia_fase_b * DR / S_base_mono; % pu

DA_fase_c = EA_dia_fase_c * DA / S_base_mono; % pu
DR_fase_c = ER_dia_fase_c * DR / S_base_mono; % pu

% Potencias aparentes de cada fase por intevalo de tempo
S_fase_a = DA_fase_a + 1j * DR_fase_a; % pu
S_fase_b = DA_fase_b + 1j * DR_fase_b; % pu
S_fase_c = DA_fase_c + 1j * DR_fase_c; % pu

% Valores base
v_ref = 1;
z_pu = (S_base_mono / V_base_fase^2) * z; % pu/km
s_a_pu = S_fase_a; % pu
s_b_pu = S_fase_b; % pu
s_c_pu = S_fase_c; % pu
v_a_e_pu = [1, 1, 1, 1]; % pu
v_b_e_pu = [-0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2]; % pu
v_c_e_pu = [-0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2]; % pu
v_n_e_pu = [0, 0, 0, 0]; % pu

% Calculos
v_a_s_pu = [1, 1, 1, 1]; % pu
v_b_s_pu = [-0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2]; % pu
v_c_s_pu = [-0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2]; % pu
v_n_s_pu = [0, 0, 0, 0]; % pu

iter = 1; % iteracao
while (erro_a(iter, :) > tolerancia) | (erro_b(iter, :) > tolerancia) | (erro_c(iter, :) > tolerancia) | (erro_n(iter, :) > tolerancia)

    i_s_a_carga_zip(iter, :) = (kp .* (conj(s_a_pu) ./ conj(v_a_s_pu(iter, :) - v_n_s_pu(iter, :)))) + (ki .* (conj(s_a_pu) ./ conj(v_ref)) .* ((v_a_s_pu(iter, :) - v_n_s_pu(iter, :)) ./ abs((v_a_s_pu(iter, :) - v_n_s_pu(iter, :))))) + (kz .* (conj(s_a_pu)./((v_ref)^2)) .* (v_a_s_pu(iter, :) - v_n_s_pu(iter, :)));
    i_s_b_carga_zip(iter, :) = (kp .* (conj(s_b_pu) ./ conj(v_b_s_pu(iter, :) - v_n_s_pu(iter, :)))) + (ki .* (conj(s_b_pu) ./ conj(v_ref)) .* ((v_b_s_pu(iter, :) - v_n_s_pu(iter, :)) ./ abs((v_b_s_pu(iter, :) - v_n_s_pu(iter, :))))) + (kz .* (conj(s_b_pu)./((v_ref)^2)) .* (v_b_s_pu(iter, :) - v_n_s_pu(iter, :)));
    i_s_c_carga_zip(iter, :) = (kp .* (conj(s_c_pu) ./ conj(v_c_s_pu(iter, :) - v_n_s_pu(iter, :)))) + (ki .* (conj(s_c_pu) ./ conj(v_ref)) .* ((v_c_s_pu(iter, :) - v_n_s_pu(iter, :)) ./ abs((v_c_s_pu(iter, :) - v_n_s_pu(iter, :))))) + (kz .* (conj(s_c_pu)./((v_ref)^2)) .* (v_c_s_pu(iter, :) - v_n_s_pu(iter, :)));
    i_s_n(iter, :) = - (i_s_a_carga_zip(iter, :) + i_s_b_carga_zip(iter, :) + i_s_c_carga_zip(iter, :));
    
    delta_v_a(iter, :) = d * z_pu * i_s_a_carga_zip(iter, :);
    delta_v_b(iter, :) = d * z_pu * i_s_b_carga_zip(iter, :);
    delta_v_c(iter, :) = d * z_pu * i_s_c_carga_zip(iter, :);
    delta_v_n(iter, :) = d * z_pu * i_s_n(iter, :);

    v_a_s_pu(iter + 1, :) = v_a_e_pu - delta_v_a(iter, :);
    v_b_s_pu(iter + 1, :) = v_b_e_pu - delta_v_b(iter, :);
    v_c_s_pu(iter + 1, :) = v_c_e_pu - delta_v_c(iter, :);
    v_n_s_pu(iter + 1, :) = v_n_e_pu - delta_v_n(iter, :);
    
    erro_a(iter + 1, :) = abs(v_a_s_pu(iter + 1, :) - v_a_s_pu(iter, :));
    erro_b(iter + 1, :) = abs(v_b_s_pu(iter + 1, :) - v_b_s_pu(iter, :));
    erro_c(iter + 1, :) = abs(v_c_s_pu(iter + 1, :) - v_c_s_pu(iter, :));
    erro_n(iter + 1, :) = abs(v_n_s_pu(iter + 1, :) - v_n_s_pu(iter, :));

    iter = iter + 1;
    
end

% Respostas
% Forma Polar
correntes_fase_a_Q1 = [round(abs(i_s_a_carga_zip), 4), round(rad2deg(angle(i_s_a_carga_zip)), 2)]
correntes_fase_b_Q1 = [round(abs(i_s_b_carga_zip), 4), round(rad2deg(angle(i_s_b_carga_zip)), 2)]
correntes_fase_c_Q1 = [round(abs(i_s_c_carga_zip), 4), round(rad2deg(angle(i_s_c_carga_zip)), 2)]
correntes_n_Q1 = [round(abs(i_s_n), 4), round(rad2deg(angle(i_s_n)), 2)]


tensoes_fase_a_Q1 = [round(abs(v_a_s_pu - v_n_s_pu), 4), round(rad2deg(angle(v_a_s_pu - v_n_s_pu)), 2)]
tensoes_fase_b_Q1 = [round(abs(v_b_s_pu - v_n_s_pu), 4), round(rad2deg(angle(v_b_s_pu - v_n_s_pu)), 2)]
tensoes_fase_c_Q1 = [round(abs(v_c_s_pu - v_n_s_pu), 4), round(rad2deg(angle(v_c_s_pu - v_n_s_pu)), 2)]
tensoes_n_Q1 = [round(abs(v_n_s_pu), 4), round(rad2deg(angle(v_n_s_pu)), 2)]

% Quest�o 2
erro_a_Q2 = [1, 1, 1, 1];
erro_b_Q2 = [1, 1, 1, 1];
erro_c_Q2 = [1, 1, 1, 1];
erro_n_Q2 = [1, 1, 1, 1];

iter_Q2 = 1; % iteracao

v_a_s_pu_Q2 = [1, 1, 1, 1]; % pu
v_b_s_pu_Q2 = [-0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2]; % pu
v_c_s_pu_Q2 = [-0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2]; % pu
v_n_s_pu_Q2 = [0, 0, 0, 0]; % pu


while (erro_a_Q2(iter_Q2, :) > tolerancia) | (erro_b_Q2(iter_Q2, :) > tolerancia) | (erro_c_Q2(iter_Q2, :) > tolerancia) | (erro_n(iter_Q2, :) > tolerancia)

    i_s_a_carga_zip_Q2(iter_Q2, :) = (kp .* (conj(s_a_pu) ./ conj(v_a_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)))) + (ki .* (conj(s_a_pu) ./ conj(v_ref)) .* ((v_a_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)) ./ abs((v_a_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :))))) + (kz .* (conj(s_a_pu)./((v_ref)^2)) .* (v_a_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)));
    i_s_b_carga_zip_Q2(iter_Q2, :) = (kp .* (conj(s_b_pu) ./ conj(v_b_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)))) + (ki .* (conj(s_b_pu) ./ conj(v_ref)) .* ((v_b_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)) ./ abs((v_b_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :))))) + (kz .* (conj(s_b_pu)./((v_ref)^2)) .* (v_b_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)));
    i_s_c_carga_zip_Q2(iter_Q2, :) = (kp .* (conj(s_c_pu) ./ conj(v_c_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)))) + (ki .* (conj(s_c_pu) ./ conj(v_ref)) .* ((v_c_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)) ./ abs((v_c_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :))))) + (kz .* (conj(s_c_pu)./((v_ref)^2)) .* (v_c_s_pu_Q2(iter_Q2, :) - v_n_s_pu_Q2(iter_Q2, :)));
    i_s_n_Q2(iter_Q2, :) = - (i_s_a_carga_zip_Q2(iter_Q2, :) + i_s_b_carga_zip_Q2(iter_Q2, :) + i_s_c_carga_zip_Q2(iter_Q2, :));
    
    delta_v_a_Q2(iter_Q2, :) = (d * z_pu * i_s_a_carga_zip_Q2(iter_Q2, :))/2;
    delta_v_b_Q2(iter_Q2, :) = (d * z_pu * i_s_b_carga_zip_Q2(iter_Q2, :))/2;
    delta_v_c_Q2(iter_Q2, :) = (d * z_pu * i_s_c_carga_zip_Q2(iter_Q2, :))/2;
    delta_v_n_Q2(iter_Q2, :) = (d * z_pu * i_s_n_Q2(iter_Q2, :))/2;

    v_a_s_pu_Q2(iter_Q2 + 1, :) = v_a_e_pu - delta_v_a_Q2(iter_Q2, :);
    v_b_s_pu_Q2(iter_Q2 + 1, :) = v_b_e_pu - delta_v_b_Q2(iter_Q2, :);
    v_c_s_pu_Q2(iter_Q2 + 1, :) = v_c_e_pu - delta_v_c_Q2(iter_Q2, :);
    v_n_s_pu_Q2(iter_Q2 + 1, :) = v_n_e_pu - delta_v_n_Q2(iter_Q2, :);
    
    erro_a_Q2(iter_Q2 + 1, :) = abs(v_a_s_pu_Q2(iter_Q2 + 1, :) - v_a_s_pu_Q2(iter_Q2, :));
    erro_b_Q2(iter_Q2 + 1, :) = abs(v_b_s_pu_Q2(iter_Q2 + 1, :) - v_b_s_pu_Q2(iter_Q2, :));
    erro_c_Q2(iter_Q2 + 1, :) = abs(v_c_s_pu_Q2(iter_Q2 + 1, :) - v_c_s_pu_Q2(iter_Q2, :));
    erro_n(iter_Q2 + 1, :) = abs(v_n_s_pu_Q2(iter_Q2 + 1, :) - v_n_s_pu_Q2(iter_Q2, :));

    iter_Q2 = iter_Q2 + 1;
    
end

% Respostas
% Forma Polar
correntes_fase_a_Q2 = [round(abs(i_s_a_carga_zip_Q2), 4), round(rad2deg(angle(i_s_a_carga_zip_Q2)), 2)]
correntes_fase_b_Q2 = [round(abs(i_s_b_carga_zip_Q2), 4), round(rad2deg(angle(i_s_b_carga_zip_Q2)), 2)]
correntes_fase_c_Q2 = [round(abs(i_s_c_carga_zip_Q2), 4), round(rad2deg(angle(i_s_c_carga_zip_Q2)), 2)]
correntes_n_Q2 = [round(abs(i_s_n_Q2), 4), round(rad2deg(angle(i_s_n_Q2)), 2)]


tensoes_fase_a_Q2 = [round(abs(v_a_s_pu_Q2 - v_n_s_pu_Q2), 4), round(rad2deg(angle(v_a_s_pu_Q2 - v_n_s_pu_Q2)), 2)]
tensoes_fase_b_Q2 = [round(abs(v_b_s_pu_Q2 - v_n_s_pu_Q2), 4), round(rad2deg(angle(v_b_s_pu_Q2 - v_n_s_pu_Q2)), 2)]
tensoes_fase_c_Q2 = [round(abs(v_c_s_pu_Q2 - v_n_s_pu_Q2), 4), round(rad2deg(angle(v_c_s_pu_Q2 - v_n_s_pu_Q2)), 2)]
tensoes_n_Q2 = [round(abs(v_n_s_pu_Q2), 4), round(rad2deg(angle(v_n_s_pu_Q2)), 2)]

