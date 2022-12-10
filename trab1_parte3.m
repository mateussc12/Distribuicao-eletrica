clc; clear;
format short
% Parametrizacao
z1 = 0.2053 + 1j*0.3753; % Ohm/km
z2 = 2 * z1; % Ohm/km
tolerancia = 0.000001;
V_se = 1; % pu

% Figura 1
d12 = 0.5; % km
d24 = 2; % km
d35 = 0.75; % km


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

% Valores base
for i = 1:length(DA)
S2_S_fases_abc(i, :) = real(S2_EA_ER) * DA(i) / S_base_mono + 1j * imag(S2_EA_ER) * DR(i) / S_base_mono; % pu
S4_S_fases_abc(i, :) = real(S4_EA_ER) * DA(i) / S_base_mono + 1j * imag(S4_EA_ER) * DR(i) / S_base_mono; % pu
S5_S_fases_abc(i, :) = real(S5_EA_ER) * DA(i) / S_base_mono + 1j * imag(S5_EA_ER) * DR(i) / S_base_mono; % pu
S6_S_fases_abc(i, :) = real(S6_EA_ER) * DA(i) / S_base_mono + 1j * imag(S6_EA_ER) * DR(i) / S_base_mono; % pu
end

z1_pu = (S_base_mono / V_base_fase^2) * z1; % pu/km
z2_pu = (S_base_mono / V_base_fase^2) * z2; % pu/km

% Questao 1

v_a_e_0 = [1, 1, 1, 1]; % pu
v_b_e_0 = [-0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2]; % pu
v_c_e_0 = [-0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2]; % pu
v_n_e_0 = [0, 0, 0, 0]; % pu

[v_a_e_2, v_b_e_2, v_c_e_2, v_n_e_2, i_s_a_2, i_s_b_2, i_s_c_2, i_s_n_2] = calcula_tensao_trifasico(S2_S_fases_abc, zip_S2, tolerancia, z1_pu, v_a_e_0, v_b_e_0, v_c_e_0, v_n_e_0, d12);
disp('Resultados:')
v_an_e_2 = abs(v_a_e_2 - v_n_e_2)
v_bn_e_2 = abs(v_b_e_2 - v_n_e_2)
v_cn_e_2 = abs(v_c_e_2 - v_n_e_2)
disp('*****************************')
[v_a_e_4, v_b_e_4, v_c_e_4, v_n_e_4, i_s_a_4, i_s_b_4, i_s_c_4, i_s_n_4] = calcula_tensao_trifasico(S2_S_fases_abc, zip_S4, tolerancia, z2_pu, v_a_e_2, v_b_e_2, v_c_e_2, v_n_e_2, d24);
%[v_a_e_5, v_b_e_5, v_c_e_5, v_n_e_5, i_s_a_5, i_s_b_5, i_s_c_5, i_s_n_5] = calcula_tensao_trifasico(S2_S_fases_abc, zip_S5, tolerancia, z2_pu, v_a_e_, v_b_e_, v_c_e_, v_n_e_, d12);

% v_a_e_2
% v_b_e_2
% v_c_e_2
% v_n_e_2

% i_s_a_2
% i_s_b_2
% i_s_c_2
% i_s_n_2
