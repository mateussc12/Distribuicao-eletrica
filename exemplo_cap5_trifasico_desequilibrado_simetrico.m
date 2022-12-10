clc; clear;
% Parametrizacao
v_base_fase = 13.8 / sqrt(3); % V
s_base = 100; % VA
z = 0.2053 + 1j*0.3753; % Ohm/km
d = 5; % km
s_a = 1 * cosd(23.074) + 1j * 1 * sind(23.074); % VA
s_b = 2 * cosd(23.074) + 1j * 2 * sind(23.074); % VA
s_c = 3 * cosd(23.074) + 1j * 3 * sind(23.074); % VA
v_a = 13.8 / sqrt(3) * cosd(0) + 1j * 13.8 / sqrt(3) * sind(0); % V
v_b = 13.8 / sqrt(3) * cosd(-120) + 1j * 13.8 / sqrt(3) * sind(-120); % V
v_c = 13.8 / sqrt(3) * cosd(120) + 1j * 13.8 / sqrt(3) * sind(120); % V
tolerancia = 0.0005;
erro_a = 1;
erro_b = 1;
erro_c = 1;
erro_n = 1;
kp = 0;
ki = 0;
kz = 1;

if kp + ki + kz > 1
    disp('kp + ki + kz devem ser <= 1')
    return
end

% Valores base
v_ref = 1;
z_pu = (s_base / v_base_fase^2) * z; % pu/km
s_a_pu = s_a / s_base; % pu
s_b_pu = s_b / s_base; % pu
s_c_pu = s_c / s_base; % pu
v_a_e_pu = v_a / v_base_fase; % pu
v_b_e_pu = v_b / v_base_fase; % pu 
v_c_e_pu = v_c / v_base_fase; % pu
v_n_e_pu = 0 + j * 0; % pu



% Calculos
v_a_s_pu = cosd(0) + 1j * sind(0); % pu
v_b_s_pu = cosd(-120) + 1j * sind(-120); % pu
v_c_s_pu = cosd(120) + 1j * sind(120); % pu
v_n_s_pu = 0 + j * 0; % pu
iter = 1; % iteracao

while (erro_a(iter) > tolerancia) || (erro_b(iter) > tolerancia) || (erro_c(iter) > tolerancia) || (erro_n(iter) > tolerancia)

    i_s_a_carga_zip(iter) = (kp * (conj(s_a_pu) / conj(v_a_s_pu(iter) - v_n_s_pu(iter)))) + (ki * (conj(s_a_pu) / conj(v_ref)) * ((v_a_s_pu(iter) - v_n_s_pu(iter)) / abs((v_a_s_pu(iter) - v_n_s_pu(iter))))) + (kz * (conj(s_a_pu)/((v_ref)^2)) * (v_a_s_pu(iter) - v_n_s_pu(iter)));
    i_s_b_carga_zip(iter) = (kp * (conj(s_b_pu) / conj(v_b_s_pu(iter) - v_n_s_pu(iter)))) + (ki * (conj(s_b_pu) / conj(v_ref)) * ((v_b_s_pu(iter) - v_n_s_pu(iter)) / abs((v_b_s_pu(iter) - v_n_s_pu(iter))))) + (kz * (conj(s_b_pu)/((v_ref)^2)) * (v_b_s_pu(iter) - v_n_s_pu(iter)));
    i_s_c_carga_zip(iter) = (kp * (conj(s_c_pu) / conj(v_c_s_pu(iter) - v_n_s_pu(iter)))) + (ki * (conj(s_c_pu) / conj(v_ref)) * ((v_c_s_pu(iter) - v_n_s_pu(iter)) / abs((v_c_s_pu(iter) - v_n_s_pu(iter))))) + (kz * (conj(s_c_pu)/((v_ref)^2)) * (v_c_s_pu(iter) - v_n_s_pu(iter)));
    i_s_n(iter) = - (i_s_a_carga_zip(iter) + i_s_b_carga_zip(iter) + i_s_c_carga_zip(iter));
    
    delta_v_a(iter) = d * z_pu * i_s_a_carga_zip(iter);
    delta_v_b(iter) = d * z_pu * i_s_b_carga_zip(iter);
    delta_v_c(iter) = d * z_pu * i_s_c_carga_zip(iter);
    delta_v_n(iter) = d * z_pu * i_s_n(iter);

    v_a_s_pu(iter + 1) = v_a_e_pu - delta_v_a(iter);
    v_b_s_pu(iter + 1) = v_b_e_pu - delta_v_b(iter);
    v_c_s_pu(iter + 1) = v_c_e_pu - delta_v_c(iter);
    v_n_s_pu(iter + 1) = v_n_e_pu - delta_v_n(iter);
    
    erro_a(iter + 1) = abs(v_a_s_pu(iter + 1) - v_a_s_pu(iter));
    erro_b(iter + 1) = abs(v_b_s_pu(iter + 1) - v_b_s_pu(iter));
    erro_c(iter + 1) = abs(v_c_s_pu(iter + 1) - v_c_s_pu(iter));
    erro_n(iter + 1) = abs(v_n_s_pu(iter + 1) - v_n_s_pu(iter));

    iter = iter + 1;
    
end

v_a_s = v_a_s_pu * v_base_fase;
abs(v_a_s(end))
rad2deg(angle(v_a_s(end)))
v_b_s = v_b_s_pu * v_base_fase;
abs(v_b_s(end))
rad2deg(angle(v_b_s(end)))
v_c_s = v_c_s_pu * v_base_fase;
abs(v_c_s(end))
rad2deg(angle(v_c_s(end)))
v_n_s = v_n_s_pu * v_base_fase;
abs(v_n_s(end))
rad2deg(angle(v_n_s(end)))
