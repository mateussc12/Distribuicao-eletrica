clc; clear; close all

V_SE = [1, 1, 1, 1, 1, 1, 1, 1.025, 1, 1, 1, 1.05, 1, 1, 1.025, 1, 1, 1, 1.025, 1.025, 1, 1, 1.025, 1.05, 1, 1, 1.05, 1, 1, 1, 1.05, 1.025, 1, 1, 1.05, 1.05, 1, 1.025, 1, 1, 1, 1.025, 1, 1.025, 1, 1.025, 1, 1.05, 1, 1.025, 1.025, 1, 1, 1.025, 1.025, 1.025, 1, 1.025, 1.025, 1.05, 1, 1.025, 1.05, 1, 1, 1.025, 1.05, 1.025, 1, 1.025, 1.05, 1.05, 1, 1.05, 1, 1, 1, 1.05, 1, 1.025, 1, 1.05, 1, 1.05, 1, 1.05, 1.025, 1, 1, 1.05, 1.025, 1.025, 1, 1.05, 1.025, 1.05, 1, 1.05, 1.05, 1, 1, 1.05, 1.05, 1.025, 1, 1.05, 1.05, 1.05, 1.025, 1, 1, 1, 1.025, 1, 1, 1.025, 1.025, 1, 1, 1.05, 1.025, 1, 1.025, 1, 1.025, 1, 1.025, 1.025, 1.025, 1, 1.025, 1.05, 1.025, 1, 1.05, 1, 1.025, 1, 1.05, 1.025, 1.025, 1, 1.05, 1.05, 1.025, 1.025, 1, 1, 1.025, 1.025, 1, 1.025, 1.025, 1.025, 1, 1.05, 1.025, 1.025, 1.025, 1, 1.025, 1.025, 1.025, 1.025, 1.025, 1.025, 1.025, 1.05, 1.025, 1.025, 1.05, 1, 1.025, 1.025, 1.05, 1.025, 1.025, 1.025, 1.05, 1.05, 1.025, 1.05, 1, 1, 1.025, 1.05, 1, 1.025, 1.025, 1.05, 1, 1.05, 1.025, 1.05, 1.025, 1, 1.025, 1.05, 1.025, 1.025, 1.025, 1.05, 1.025, 1.05, 1.025, 1.05, 1.05, 1, 1.025, 1.05, 1.05, 1.025, 1.025, 1.05, 1.05, 1.05, 1.05, 1, 1, 1, 1.05, 1, 1, 1.025, 1.05, 1, 1, 1.05, 1.05, 1, 1.025, 1, 1.05, 1, 1.025, 1.025, 1.05, 1, 1.025, 1.05, 1.05, 1, 1.05, 1, 1.05, 1, 1.05, 1.025, 1.05, 1, 1.05, 1.05, 1.05, 1.025, 1, 1, 1.05, 1.025, 1, 1.025, 1.05, 1.025, 1, 1.05, 1.05, 1.025, 1.025, 1, 1.05, 1.025, 1.025, 1.025, 1.05, 1.025, 1.025, 1.05, 1.05, 1.025, 1.05, 1, 1.05, 1.025, 1.05, 1.025, 1.05, 1.025, 1.05, 1.05, 1.05, 1.05, 1, 1, 1.05, 1.05, 1, 1.025, 1.05, 1.05, 1, 1.05, 1.05, 1.05, 1.025, 1, 1.05, 1.05, 1.025, 1.025, 1.05, 1.05, 1.025, 1.05, 1.05, 1.05, 1.05, 1, 1.05, 1.05, 1.05, 1.025, 1.05, 1.05, 1.05, 1.05];
custo_otimo = 999999; % reais
condition = 1;
i = 1;
while condition
   
	if i+3 >= length(V_SE)
		condition = 0;
	end

	for banco_cap = 1:5
		[table_perdas_tecnicas, v2, v4, v5, v6, table_custo_operacional, custo_operacional_total] = calc_fluxo(V_SE(i:(i+3)), banco_cap);

		if custo_operacional_total <= custo_otimo
			custo_otimo = custo_operacional_total;
			v_se_otimo = V_SE(i:(i+3));
	        banco_otimo = banco_cap;

			table_perdas_tecnicas_otima = table_perdas_tecnicas;
			v2an_otima = v2(1:8);
			v2bn_otima = v2(9:16);
			v2cn_otima = v2(17:24);
			v4an_otima = v4(1:8);
			v4bn_otima = v4(9:16);
			v4cn_otima = v4(17:24);
			v5an_otima = v5(1:8);
			v5bn_otima = v5(9:16);
			v5cn_otima = v5(17:24);
			v6an_otima = v6(1:8);
			v6bn_otima = v6(9:16);
			v6cn_otima = v6(17:24);
            table_custo_operacional_otima = table_custo_operacional;
            
		end
	end
i = i + 4;
end

% Tabelas
table_perdas_tecnicas_otima
table_custo_operacional_otima
v_se_otimo
custo_otimo
if banco_otimo == 1
	banco = "Sem banco de capacitores"
elseif banco_otimo == 2
	banco = "Banco de capacitores na barra 2"
elseif banco_otimo == 3
	banco = "Banco de capacitores na barra 4"
elseif banco_otimo == 4
	banco = "Banco de capacitores na barra 5"
elseif banco_otimo == 5
	banco = "Banco de capacitores na barra 6"
end



% Graficos
t = [0, 0.25, 0.5, 0.75, 1];

figure
plot_v2an = [v2an_otima(1), v2an_otima(2), v2an_otima(3), v2an_otima(4), v2an_otima(1)];
plot_v2bn = [v2bn_otima(1), v2bn_otima(2), v2bn_otima(3), v2bn_otima(4), v2bn_otima(1)];
plot_v2cn = [v2cn_otima(1), v2cn_otima(2), v2cn_otima(3), v2cn_otima(4), v2cn_otima(1)];
hold all
stairs(t, plot_v2an, 'linewidth', 2, 'color', 'g')
stairs(t, plot_v2bn, 'linewidth', 2, 'color', 'black')
stairs(t, plot_v2cn, 'linewidth', 2, 'color', 'm')
plot(t, [0.9 0.9 0.9 0.9 0.9], '--', 'color', 'r', 'linewidth', 2)
plot(t, [0.93 0.93 0.93 0.93 0.93], '--', 'color', 'b', 'linewidth', 2)
plot(t, [1.05 1.05 1.05 1.05 1.05], '--', 'color', 'r', 'linewidth', 2)
title("V2an")
ylabel("Amplitude (pu)")
xlabel("Tempo (%)")
legend(["Fase A", "Fase B", "Fase C"])

figure
plot_v4an = [v4an_otima(1), v4an_otima(2), v4an_otima(3), v4an_otima(4), v4an_otima(1)];
plot_v4bn = [v4bn_otima(1), v4bn_otima(2), v4bn_otima(3), v4bn_otima(4), v4bn_otima(1)];
plot_v4cn = [v4cn_otima(1), v4cn_otima(2), v4cn_otima(3), v4cn_otima(4), v4cn_otima(1)];
hold all
stairs(t, plot_v4an, 'linewidth', 2, 'color', 'g')
stairs(t, plot_v4bn, 'linewidth', 2, 'color', 'black')
stairs(t, plot_v4cn, 'linewidth', 2, 'color', 'm')
plot(t, [0.9 0.9 0.9 0.9 0.9], '--', 'color', 'r', 'linewidth', 2)
plot(t, [0.93 0.93 0.93 0.93 0.93], '--', 'color', 'b', 'linewidth', 2)
plot(t, [1.05 1.05 1.05 1.05 1.05], '--', 'color', 'r', 'linewidth', 2)
title("V4an")
ylabel("Amplitude (pu)")
xlabel("Tempo (%)")
legend(["Fase A", "Fase B", "Fase C"])

figure
plot_v5an = [v5an_otima(1), v5an_otima(2), v5an_otima(3), v5an_otima(4), v5an_otima(1)];
plot_v5bn = [v5bn_otima(1), v5bn_otima(2), v5bn_otima(3), v5bn_otima(4), v5bn_otima(1)];
plot_v5cn = [v5cn_otima(1), v5cn_otima(2), v5cn_otima(3), v5cn_otima(4), v5cn_otima(1)];
hold all
stairs(t, plot_v5an, 'linewidth', 2, 'color', 'g')
stairs(t, plot_v5bn, 'linewidth', 2, 'color', 'black')
stairs(t, plot_v5cn, 'linewidth', 2, 'color', 'm')
plot(t, [0.9 0.9 0.9 0.9 0.9], '--', 'color', 'r', 'linewidth', 2)
plot(t, [0.93 0.93 0.93 0.93 0.93], '--', 'color', 'b', 'linewidth', 2)
plot(t, [1.05 1.05 1.05 1.05 1.05], '--', 'color', 'r', 'linewidth', 2)
title("V5an")
ylabel("Amplitude (pu)")
xlabel("Tempo (%)")
legend(["Fase A", "Fase B", "Fase C"])

figure
plot_v6an = [v6an_otima(1), v6an_otima(2), v6an_otima(3), v6an_otima(4), v6an_otima(1)];
plot_v6bn = [v6bn_otima(1), v6bn_otima(2), v6bn_otima(3), v6bn_otima(4), v6bn_otima(1)];
plot_v6cn = [v6cn_otima(1), v6cn_otima(2), v6cn_otima(3), v6cn_otima(4), v6cn_otima(1)];
hold all
stairs(t, plot_v6an, 'linewidth', 2, 'color', 'g')
stairs(t, plot_v6bn, 'linewidth', 2, 'color', 'black')
stairs(t, plot_v6cn, 'linewidth', 2, 'color', 'm')
plot(t, [0.9 0.9 0.9 0.9 0.9], '--', 'color', 'r', 'linewidth', 2)
plot(t, [0.93 0.93 0.93 0.93 0.93], '--', 'color', 'b', 'linewidth', 2)
plot(t, [1.05 1.05 1.05 1.05 1.05], '--', 'color', 'r', 'linewidth', 2)
title("V6an")
ylabel("Amplitude (pu)")
xlabel("Tempo (%)")
legend(["Fase A", "Fase B", "Fase C"])
