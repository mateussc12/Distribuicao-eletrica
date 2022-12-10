function [v_a_e, v_b_e, v_c_e, v_n_e, i_s_a, i_s_b, i_s_c, i_s_n] = calcula_tensao_trifasico(S_fases_abc, modelo_zip, tolerancia, z_pu, v_a_e_pu, v_b_e_pu, v_c_e_pu, v_n_e_pu, d)

	v_a_s_pu = [1, 1, 1, 1]; % pu
	v_b_s_pu = [-0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2, -0.5 - 1j * sqrt(3) / 2]; % pu
	v_c_s_pu = [-0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2, -0.5 + 1j * sqrt(3) / 2]; % pu
	v_n_s_pu = [0, 0, 0, 0]; % pu

	v_ref = 1; % pu

	for i = 1:3
		S_fase_a = S_fases_abc(:, 1);
		S_fase_b = S_fases_abc(:, 2);
		S_fase_c = S_fases_abc(:, 3);
	end

	S_fase_a = transpose(S_fase_a);
	S_fase_b = transpose(S_fase_b);
	S_fase_c = transpose(S_fase_c);
    
	erro_a = [1, 1, 1, 1];
	erro_b = [1, 1, 1, 1];
	erro_c = [1, 1, 1, 1];
	erro_n = [1, 1, 1, 1];

	iter = 1;
while (erro_a(iter, :) > tolerancia) | (erro_b(iter, :) > tolerancia) | (erro_c(iter, :) > tolerancia) | (erro_n(iter, :) > tolerancia)
    iter
    v_a_s_pu(iter, :)

    i_s_a_carga_zip(iter, :) = (modelo_zip(3) .* (conj(S_fase_a) ./ conj(v_a_s_pu(iter, :) - v_n_s_pu(iter, :)))) + (modelo_zip(2) .* (conj(S_fase_a) ./ conj(v_ref)) .* ((v_a_s_pu(iter, :) - v_n_s_pu(iter, :)) ./ abs((v_a_s_pu(iter, :) - v_n_s_pu(iter, :))))) + (modelo_zip(1) .* (conj(S_fase_a)./((v_ref)^2)) .* (v_a_s_pu(iter, :) - v_n_s_pu(iter, :)));
    i_s_a_carga_zip(iter, :)
    i_s_b_carga_zip(iter, :) = (modelo_zip(3) .* (conj(S_fase_b) ./ conj(v_b_s_pu(iter, :) - v_n_s_pu(iter, :)))) + (modelo_zip(2) .* (conj(S_fase_b) ./ conj(v_ref)) .* ((v_b_s_pu(iter, :) - v_n_s_pu(iter, :)) ./ abs((v_b_s_pu(iter, :) - v_n_s_pu(iter, :))))) + (modelo_zip(1) .* (conj(S_fase_b)./((v_ref)^2)) .* (v_b_s_pu(iter, :) - v_n_s_pu(iter, :)));
    i_s_c_carga_zip(iter, :) = (modelo_zip(3) .* (conj(S_fase_c) ./ conj(v_c_s_pu(iter, :) - v_n_s_pu(iter, :)))) + (modelo_zip(2) .* (conj(S_fase_c) ./ conj(v_ref)) .* ((v_c_s_pu(iter, :) - v_n_s_pu(iter, :)) ./ abs((v_c_s_pu(iter, :) - v_n_s_pu(iter, :))))) + (modelo_zip(1) .* (conj(S_fase_c)./((v_ref)^2)) .* (v_c_s_pu(iter, :) - v_n_s_pu(iter, :)));
    i_s_n_carga_zip(iter, :) = - (i_s_a_carga_zip(iter, :) + i_s_b_carga_zip(iter, :) + i_s_c_carga_zip(iter, :));
    
    delta_v_a(iter, :) = (d * z_pu * i_s_a_carga_zip(iter, :))/2;
    delta_v_b(iter, :) = (d * z_pu * i_s_b_carga_zip(iter, :))/2;
    delta_v_c(iter, :) = (d * z_pu * i_s_c_carga_zip(iter, :))/2;
    delta_v_n(iter, :) = (d * z_pu * i_s_n_carga_zip(iter, :))/2;
    delta_v_a(iter, :)

    v_a_s_pu(iter + 1, :) = v_a_e_pu - delta_v_a(iter, :);
    v_a_s_pu(iter + 1, :)
    disp('####################################')
    v_b_s_pu(iter + 1, :) = v_b_e_pu - delta_v_b(iter, :);
    v_c_s_pu(iter + 1, :) = v_c_e_pu - delta_v_c(iter, :);
    v_n_s_pu(iter + 1, :) = v_n_e_pu - delta_v_n(iter, :);

    erro_a(iter + 1, :) = abs(v_a_s_pu(iter + 1, :) - v_a_s_pu(iter, :));
    erro_b(iter + 1, :) = abs(v_b_s_pu(iter + 1, :) - v_b_s_pu(iter, :));
    erro_c(iter + 1, :) = abs(v_c_s_pu(iter + 1, :) - v_c_s_pu(iter, :));
    erro_n(iter + 1, :) = abs(v_n_s_pu(iter + 1, :) - v_n_s_pu(iter, :));

    iter = iter + 1;

end
    
    v_a_e = v_a_s_pu(end, :);
    v_b_e = v_b_s_pu(end, :); 
    v_c_e = v_c_s_pu(end, :);
    v_n_e = v_n_s_pu(end, :);
    i_s_a = i_s_a_carga_zip(end, :);
    i_s_b = i_s_b_carga_zip(end, :);
    i_s_c = i_s_c_carga_zip(end, :);
    i_s_n = i_s_n_carga_zip(end, :);

end