function [S_fase_a, S_fase_b, S_fase_c] = ordem_potencia(S_fases_abc)

	for i = 1:3
		S_fase_a = S_fases_abc(:, 1);
		S_fase_b = S_fases_abc(:, 2);
		S_fase_c = S_fases_abc(:, 3);
	end

	S_fase_a = transpose(S_fase_a);
	S_fase_b = transpose(S_fase_b);
	S_fase_c = transpose(S_fase_c);

end