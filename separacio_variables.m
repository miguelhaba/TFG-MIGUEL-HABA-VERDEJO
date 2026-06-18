clear; clc; close all;

% Constants
alpha = 7.2973525693e-3;
c = 299792458e15; % fm/s
h = 6.62607015e-34 * 1.602176634e-19 / (2*pi); % eV·s
hc = h * c; % eV·fm
me = 0.51099895000e6; % eV
a0 = hc / (me * alpha); % fm
E1s = -alpha * hc / (2 * a0); % eV

mp = 938.27208816e6; % eV
mu = mp / 2; % eV 
K = hc^2 / (2 * mu); % eV·fm²

%% CÀLCUL ENERGIA
dosRa0 = linspace(0.3, 10, 10000);
R_fm = dosRa0/2 * a0;  % R en fm

% Càlcul LCAO
[E_LCAO, S, H11, H12] = calcul_LCAO(dosRa0, alpha, hc, a0, E1s);

% Càlcul pel mètode de separació de variables (N=4)
N = 4;
[E_exacta, C_opt, info_exacte] = calcul_separacio_variables(dosRa0, R_fm, E_LCAO, alpha, hc, me, N);

fprintf('\n========== RESULTATS DEL MÈTODE DE SEPARACIÓ DE VARIABLES ==========\n');
fprintf('Distància d''equilibri: 2R = %.4f a0\n', info_exacte.dosRa0_min);
fprintf('Energia mínima exacta: E = %.7f eV\n', info_exacte.E_min);
fprintf('====================================================================\n\n');

% 2.3. Mínim LCAO
[E_min_LCAO, idx_LCAO] = min(E_LCAO);
R_min_LCAO = dosRa0(idx_LCAO);

%% GRÀFICS

% Gràfic convergència
figure('Name', 'Convergència del mètode de separació de variables', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);
hold on;

% Colors per als diferents ordres
colors = lines(4);
colors(3, :) = [0.10 0.55 0.25];

% Calcular per a N=1,2,3,4 i dibuixar
for N_plot = 1:4
    [E_plot, ~, info_plot] = calcul_separacio_variables(dosRa0, R_fm, E_LCAO, alpha, hc, me, N_plot);
    
    plot(dosRa0, E_plot, 'Color', colors(N_plot, :), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('N=%d. Mínim: (%.4fa_0, %.7f eV)', ...
         N_plot, info_plot.dosRa0_min, info_plot.E_min));
    plot(info_plot.dosRa0_min, info_plot.E_min, 'o', ...
         'MarkerEdgeColor', colors(N_plot, :), ...
         'MarkerFaceColor', 'none', 'MarkerSize', 8, ...
         'HandleVisibility', 'off');
end

% LCAO i E1s
plot(dosRa0, E_LCAO, '--k', 'LineWidth', 1.2, ...
     'DisplayName', sprintf('LCAO. Mínim: (%.4fa_0, %.4f eV)', R_min_LCAO, E_min_LCAO));
plot([0 10], [E1s E1s], 'k-.', 'DisplayName', 'Energia E_{1s}');
plot(R_min_LCAO, E_min_LCAO, 'dk', 'MarkerFaceColor', 'none', ...
     'MarkerSize', 8, 'HandleVisibility', 'off');

xlabel('Distància internuclear (2R/a_0)', 'FontSize', 12);
ylabel('Energia (eV)', 'FontSize', 12);
title('Convergència del mètode de separació de variables', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 12); 
ylim([-17, -12.5]); xlim([0, 10]);
grid off;
hold off;

%% FUNCIÓ D'ONA
[Psi_exacta, X, Z, info_ona] = funcio_ona_2D(info_exacte.dosRa0_min, info_exacte.R_min, info_exacte.C_min, info_exacte.B_min, E1s);

% Gràfic 2D de la funció d'ona exacta
figure('Position', [150, 150, 1100, 500],'Name', 'Funció d''ona - Mètode de separació de variables', 'NumberTitle', 'off');

contourf(X, Z, Psi_exacta, 20);
colorbar;
xlabel('x (a_0)', 'FontSize', 12);
ylabel('z (a_0)', 'FontSize', 12);
title(sprintf('Funció d''ona \\phi (2R = %.4f a_0)', info_exacte.dosRa0_min), 'FontSize', 14);
axis equal;
colormap(jet);

%% NUMEROV
resultats_num = numerov(dosRa0, R_fm, E_exacta, a0, K);

fprintf('\n========== ESTATS VIBRACIONALS ==========\n');
for iv = 0:2
    fprintf('Energia v=%d: %.7f eV\n', iv, resultats_num.v{iv+1}.E);
end

% Gràfica del potencial amb nivells vibracionals
figure('Name', 'Estats vibracionals', 'NumberTitle', 'off', ...
       'Position', [100, 100, 900, 600]);
hold on;

c_v0 = [0.70 0.15 0.15];
c_v1 = [0.00 0.35 0.60];
c_v2 = [0.75 0.60 0.10];

% Regions 
punts_tall0 = find(E_exacta <= resultats_num.v{1}.E);
punts_tall1 = find(E_exacta <= resultats_num.v{2}.E);
punts_tall2 = find(E_exacta <= resultats_num.v{3}.E);

fill([dosRa0(punts_tall2), fliplr(dosRa0(punts_tall2))], ...
     [E_exacta(punts_tall2), fliplr(resultats_num.v{3}.E * ones(size(punts_tall2)))], ...
     c_v2, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([dosRa0(punts_tall1), fliplr(dosRa0(punts_tall1))], ...
     [E_exacta(punts_tall1), fliplr(resultats_num.v{2}.E * ones(size(punts_tall1)))], ...
     c_v1, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([dosRa0(punts_tall0), fliplr(dosRa0(punts_tall0))], ...
     [E_exacta(punts_tall0), fliplr(resultats_num.v{1}.E * ones(size(punts_tall0)))], ...
     c_v0, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(dosRa0, E_exacta, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Potencial electrònic exacte');
plot(dosRa0, E_LCAO, 'k--', 'DisplayName', 'Aproximació LCAO');
plot([0 10], [E1s E1s], 'k-.', 'DisplayName', 'Energia E_{1s}');
plot(info_exacte.dosRa0_min, info_exacte.E_min, 'bo', ...
     'MarkerSize', 6, 'MarkerFaceColor', 'none', 'DisplayName', ...
     sprintf('Mínim exacte: (%.4fa_0, %.4f eV)', info_exacte.dosRa0_min, info_exacte.E_min));
plot(R_min_LCAO, E_min_LCAO, 'kd', 'MarkerSize', 6, 'DisplayName', ...
     sprintf('Mínim LCAO: (%.4fa_0, %.4f eV)', R_min_LCAO, E_min_LCAO));

% Línies de nivell
h_v0 = line([dosRa0(punts_tall0(1)), dosRa0(punts_tall0(end))], ...
            [resultats_num.v{1}.E, resultats_num.v{1}.E], ...
            'Color', c_v0, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('v=0: %.4f eV', resultats_num.v{1}.E));
h_v1 = line([dosRa0(punts_tall1(1)), dosRa0(punts_tall1(end))], ...
            [resultats_num.v{2}.E, resultats_num.v{2}.E], ...
            'Color', c_v1, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('v=1: %.4f eV', resultats_num.v{2}.E));
h_v2 = line([dosRa0(punts_tall2(1)), dosRa0(punts_tall2(end))], ...
            [resultats_num.v{3}.E, resultats_num.v{3}.E], ...
            'Color', c_v2, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('v=2: %.4f eV', resultats_num.v{3}.E));

xlabel('Distància internuclear (2R/a_0)', 'FontSize', 12);
ylabel('Energia (eV)', 'FontSize', 12);
title('Potencial electrònic i estats lligats del H_2^+', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 12); 
xlim([0.5, 7]); ylim([-17, -12]);
grid off;
hold off;

%% NIVELLS ROTACIONALS
fprintf('\n\n\n================ ESQUEMA D''ENERGIES ROTACIONALS (eV) ================\n');
fprintf('%-6s | %-12s | %-12s | %-12s | %-12s\n', 'v \\ L', 'L = 0', 'L = 1', 'L = 2', 'L = 3');
fprintf('-------------------------------------------------------------------------\n');
for v_idx = 1:3
    fprintf('v = %d   | %-12.7f | %-12.7f | %-12.7f | %-12.7f\n', ...
        v_idx-1, resultats_num.Energies_vL(v_idx,1), resultats_num.Energies_vL(v_idx,2), ...
        resultats_num.Energies_vL(v_idx,3), resultats_num.Energies_vL(v_idx,4));
end

% Diagrama de nivells
figure('Name', 'Diagrama de nivells rotacionals', 'NumberTitle', 'off', ...
       'Position', [100, 100, 600, 800]);
hold on;

colors_L = [0 0 0; 0 0.2 0.6; 0.6 0 0; 0 0.5 0];
x_esq = 0; x_dret = 1;

for v_idx = 1:3
    for l_idx = 1:length(resultats_num.L_valors)
        E_nivell = resultats_num.Energies_vL(v_idx, l_idx);
        plot([x_esq, x_dret], [E_nivell, E_nivell], ...
             'Color', colors_L(l_idx, :), 'LineWidth', 1.5);
    end
    text(x_esq - 0.05, resultats_num.Energies_vL(v_idx, 1), ...
         sprintf('v = %d', v_idx-1), 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

ylabel('Energia (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('Diagrama energètic del H_2^+', 'FontSize', 14);
xlim([-0.5, 1.5]);
ylim([min(resultats_num.Energies_vL(:)) - 0.2, max(resultats_num.Energies_vL(:)) + 0.2]);

h_leg = zeros(1, 4);
for l_idx = 1:4
    h_leg(l_idx) = plot(nan, nan, 'Color', colors_L(l_idx, :), 'LineWidth', 1.5);
end
legend(h_leg, {'L = 0', 'L = 1', 'L = 2', 'L = 3'}, 'Location', 'northwest', 'FontSize', 12);
ax = gca; ax.XAxis.Visible = 'off';
hold off;