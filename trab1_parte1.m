clc; clear;
% Distribuicao - Trab 1
% Mateus Souza Coelho					
% Catarina Sastre Neves

%##########################Enunciado#######################################
% 1 Tabela
EA_mes_fase_a = 1000e6; % Wh/mes
EA_mes_fase_b = 400e6; % Wh/mes
EA_mes_fase_c = 1200e6; % Wh/mes
ER_mes_fase_a = 700e6; % Wh/mes
ER_mes_fase_b = 400e6; % Wh/mes
ER_mes_fase_c = 900e6; % Wh/mes
% 2 Tabela
DA = [0.5, 0.8, 1.3, 0.7];
DR = [0.4, 0.6, 1, 0.5];
% 3 tabela
Kz = 0.5;
Ki = 0.2;
Kp = 0.3;
% Parametrizacao
S_base_tri = 10e6; % VA
S_base_mono = S_base_tri / 3; % VA
V_base_linha = 11.4e3; % V
V_base_fase = V_base_linha / sqrt(3); % V
%##########################################################################

%##########################Questao 1#######################################
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


vetor_tempo = 0:6:24;

hold on;
stairs(vetor_tempo, [DA_fase_a, DA_fase_a(end)], 'linewidth', 2);
stairs(vetor_tempo, [DA_fase_b, DA_fase_b(end)], 'linewidth', 2);
stairs(vetor_tempo, [DA_fase_c, DA_fase_c(end)], 'linewidth', 2);
title('Demandas ativas por fase')
ylabel('pu')
xlabel('Horas')
legend('Fase a', 'Fase b', 'Fase c')
xticks([0, 6, 12, 18, 24]);

figure;
hold on;
stairs(vetor_tempo, [DR_fase_a, DR_fase_a(end)], 'linewidth', 2);
stairs(vetor_tempo, [DR_fase_b, DR_fase_b(end)], 'linewidth', 2);
stairs(vetor_tempo, [DR_fase_c, DR_fase_c(end)], 'linewidth', 2);
title('Demandas reativas por fase')
ylabel('pu')
xlabel('Horas')
legend('Fase a', 'Fase b', 'Fase c')
xticks([0, 6, 12, 18, 24]);
%##########################################################################

%##########################Questao 2#######################################
% Potencias aparentes de cada fase por intevalo de tempo
S_fase_a = DA_fase_a + 1j * DR_fase_a; % pu
S_fase_b = DA_fase_b + 1j * DR_fase_b; % pu
S_fase_c = DA_fase_c + 1j * DR_fase_c; % pu

V_nominal_pu_Q2 = 1;
V_an_Q2 = 1; % 1<0 pu
V_bn_Q2 = -0.5 - 1j * sqrt(3) / 2; % 1<-120 pu
V_cn_Q2 = -0.5 + 1j * sqrt(3) / 2; % 1<120 pu

corrente_fase_a_Q2 = (Kp * conj(S_fase_a) / conj(V_an_Q2)) + (Ki * conj(S_fase_a) / V_nominal_pu_Q2 * (V_an_Q2 / abs(V_an_Q2))) + (Kz * conj(S_fase_a) / V_nominal_pu_Q2^2 * V_an_Q2); % pu
corrente_fase_b_Q2 = (Kp * conj(S_fase_b) / conj(V_bn_Q2)) + (Ki * conj(S_fase_b) / V_nominal_pu_Q2 * (V_bn_Q2 / abs(V_bn_Q2))) + (Kz * conj(S_fase_b) / V_nominal_pu_Q2^2 * V_bn_Q2); % pu
corrente_fase_c_Q2 = (Kp * conj(S_fase_c) / conj(V_cn_Q2)) + (Ki * conj(S_fase_c) / V_nominal_pu_Q2 * (V_cn_Q2 / abs(V_cn_Q2))) + (Kz * conj(S_fase_c) / V_nominal_pu_Q2^2 * V_cn_Q2); % pu
corrente_neutro_Q2 = - (corrente_fase_a_Q2 + corrente_fase_b_Q2 + corrente_fase_c_Q2); % pu
%##########################################################################

%##########################Questao 3#######################################
V_nominal_pu_Q3 = 1;
V_an_Q3 = 1.1; % 1.1<0 pu
V_bn_Q3 = -0.55 - 1j * 11 * sqrt(3) / 20; % 1.1<-120 pu
V_cn_Q3 = -0.55 + 1j * 11 * sqrt(3) / 20; % 1.1<120 pu

corrente_fase_a_Q3 = (Kp * conj(S_fase_a) / conj(V_an_Q3)) + (Ki * conj(S_fase_a) / V_nominal_pu_Q3 * (V_an_Q3 / abs(V_an_Q3))) + (Kz * conj(S_fase_a) / V_nominal_pu_Q3^2 * V_an_Q3); % pu
corrente_fase_b_Q3 = (Kp * conj(S_fase_b) / conj(V_bn_Q3)) + (Ki * conj(S_fase_b) / V_nominal_pu_Q3 * (V_bn_Q3 / abs(V_bn_Q3))) + (Kz * conj(S_fase_b) / V_nominal_pu_Q3^2 * V_bn_Q3); % pu
corrente_fase_c_Q3 = (Kp * conj(S_fase_c) / conj(V_cn_Q3)) + (Ki * conj(S_fase_c) / V_nominal_pu_Q3 * (V_cn_Q3 / abs(V_cn_Q3))) + (Kz * conj(S_fase_c) / V_nominal_pu_Q3^2 * V_cn_Q3); % pu
corrente_neutro_Q3 = - (corrente_fase_a_Q3 + corrente_fase_b_Q3 + corrente_fase_c_Q3); % pu
%##########################################################################

%##########################Questao 4#######################################
V_nominal_pu_Q4 = 1;
V_an_Q4 = 0.9; % 0.9<0 pu
V_bn_Q4 = -0.45 - 1j * 9 * sqrt(3) / 20; % 0.9<-120 pu
V_cn_Q4 = -0.45 + 1j * 9 * sqrt(3) / 20; % 0.9<120 pu

corrente_fase_a_Q4 = (Kp * conj(S_fase_a) / conj(V_an_Q4)) + (Ki * conj(S_fase_a) / V_nominal_pu_Q4 * (V_an_Q4 / abs(V_an_Q4))) + (Kz * conj(S_fase_a) / V_nominal_pu_Q4^2 * V_an_Q4); % pu
corrente_fase_b_Q4 = (Kp * conj(S_fase_b) / conj(V_bn_Q4)) + (Ki * conj(S_fase_b) / V_nominal_pu_Q4 * (V_bn_Q4 / abs(V_bn_Q4))) + (Kz * conj(S_fase_b) / V_nominal_pu_Q4^2 * V_bn_Q4); % pu
corrente_fase_c_Q4 = (Kp * conj(S_fase_c) / conj(V_cn_Q4)) + (Ki * conj(S_fase_c) / V_nominal_pu_Q4 * (V_cn_Q4 / abs(V_cn_Q4))) + (Kz * conj(S_fase_c) / V_nominal_pu_Q4^2 * V_cn_Q4); % pu
corrente_neutro_Q4 = - (corrente_fase_a_Q4 + corrente_fase_b_Q4 + corrente_fase_c_Q4); % pu
%##########################################################################

%##########################Questao 5#######################################
V_nominal_pu_Q5 = 1;
V_an_Q5 = 0.9;
V_bn_Q5 = -0.5 - 1j * sqrt(3) / 2; % 1<-120 pu
V_cn_Q5 = -0.5 + 1j * sqrt(3) / 2; % 1<120 pu

corrente_fase_a_Q5 = (Kp * conj(S_fase_a) / conj(V_an_Q5)) + (Ki * conj(S_fase_a) / V_nominal_pu_Q5 * (V_an_Q5 / abs(V_an_Q5))) + (Kz * conj(S_fase_a) / V_nominal_pu_Q5^2 * V_an_Q5); % pu
corrente_fase_b_Q5 = (Kp * conj(S_fase_b) / conj(V_bn_Q5)) + (Ki * conj(S_fase_b) / V_nominal_pu_Q5 * (V_bn_Q5 / abs(V_bn_Q5))) + (Kz * conj(S_fase_b) / V_nominal_pu_Q5^2 * V_bn_Q5); % pu
corrente_fase_c_Q5 = (Kp * conj(S_fase_c) / conj(V_cn_Q5)) + (Ki * conj(S_fase_c) / V_nominal_pu_Q5 * (V_cn_Q5 / abs(V_cn_Q5))) + (Kz * conj(S_fase_c) / V_nominal_pu_Q5^2 * V_cn_Q5); % pu
corrente_neutro_Q5 = - (corrente_fase_a_Q5 + corrente_fase_b_Q5 + corrente_fase_c_Q5); % pu
%##########################################################################

%##########################Formatacao respostas############################
% Forma Retangular
D_Q1 = [
    round(DA_fase_a, 4) + 1j * round(DR_fase_a, 4);
    round(DA_fase_b, 4) + 1j * round(DR_fase_b, 4);
    round(DA_fase_c, 4) + 1j * round(DR_fase_c, 4);
    ]
    
% Forma Polar
correntes_Q2 = [ 
    round(abs(corrente_fase_a_Q2), 4), round(rad2deg(angle(corrente_fase_a_Q2)), 2);
    round(abs(corrente_fase_b_Q2), 4), round(rad2deg(angle(corrente_fase_b_Q2)), 2);
    round(abs(corrente_fase_c_Q2), 4), round(rad2deg(angle(corrente_fase_c_Q2)), 2);
    round(abs(corrente_neutro_Q2), 4), round(rad2deg(angle(corrente_neutro_Q2)), 2)
    ]

correntes_Q3 = [ 
    round(abs(corrente_fase_a_Q3), 4), round(rad2deg(angle(corrente_fase_a_Q3)), 2);
    round(abs(corrente_fase_b_Q3), 4), round(rad2deg(angle(corrente_fase_b_Q3)), 2);
    round(abs(corrente_fase_c_Q3), 4), round(rad2deg(angle(corrente_fase_c_Q3)), 2);
    round(abs(corrente_neutro_Q3), 4), round(rad2deg(angle(corrente_neutro_Q3)), 2)
    ]

correntes_Q4 = [ 
    round(abs(corrente_fase_a_Q4), 4), round(rad2deg(angle(corrente_fase_a_Q4)), 2);
    round(abs(corrente_fase_b_Q4), 4), round(rad2deg(angle(corrente_fase_b_Q4)), 2);
    round(abs(corrente_fase_c_Q4), 4), round(rad2deg(angle(corrente_fase_c_Q4)), 2);
    round(abs(corrente_neutro_Q4), 4), round(rad2deg(angle(corrente_neutro_Q4)), 2)
    ]

correntes_Q5 = [ 
    round(abs(corrente_fase_a_Q5), 4), round(rad2deg(angle(corrente_fase_a_Q5)), 2);
    round(abs(corrente_fase_b_Q5), 4), round(rad2deg(angle(corrente_fase_b_Q5)), 2);
    round(abs(corrente_fase_c_Q5), 4), round(rad2deg(angle(corrente_fase_c_Q5)), 2);
    round(abs(corrente_neutro_Q5), 4), round(rad2deg(angle(corrente_neutro_Q5)), 2)
    ]
%##########################################################################