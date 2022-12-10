clc; clear;
% Parametrizacao
v_base = 13800; % V
s_base = 100e6; % VA
z = 0.2053 + 1j*0.3753; % Ohm/km
d = 5; % km
s = 6e6 * cosd(23.074) + 1j * 6e6 * sind(23.074); % VA
v = 13800; % V
tolerancia = 0.0005;
erro = 1;
kp = 0;
ki = 0;
kz = 1;

if kp + ki + kz > 1
    disp('kp + ki + kz devem ser <= 1')
    return
end

% Valores base
z_pu = (s_base / v_base^2) * z; % pu/km
s_pu = s / s_base; % pu
v_e_pu = v / v_base; % pu 


% Calculos
v_s_pu = 1; % pu
iter = 1; % iteracao
while erro(iter) > tolerancia

    i_s_carga_zip(iter) = kp * (conj(s_pu)/conj(v_s_pu(iter))) + ki * (conj(s_pu)/conj(v_e_pu)) * (v_s_pu(iter)/abs(v_s_pu(iter))) + kz * (conj(s_pu)/((v_e_pu)^2)) * v_s_pu(iter);

    delta_v(iter) = d * z_pu * i_s_carga_zip(iter);

    v_s_pu(iter + 1) = v_e_pu - delta_v(iter);

    erro(iter + 1) = abs(v_s_pu(iter + 1) - v_s_pu(iter));

    iter = iter + 1;

end

v_s = v_s_pu * v_base;
abs(v_s(end))
rad2deg(angle(v_s(end)))
