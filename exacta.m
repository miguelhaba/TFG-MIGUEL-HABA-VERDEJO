%% Mètode exacte i Numerov 
clear; clc; close all;

% Constants
alpha = 7.2973525693e-3; 
hc = 197.3269804e6; % eV·fm
me = 0.51099895e6; % eV
a0 = hc/(me*alpha); % fm
E1s = -alpha*hc/2/a0; % eV

dosRa0 = linspace(0.3, 10, 10000);
R = dosRa0/2 * a0;

figure(1)
hold on;
colors = lines(4);
colors(3, :) = [0.10 0.55 0.25]; 

E_tots_ordres = zeros(4, length(R));
R_tots_ordres = zeros(4, 1);
e_R = zeros(4,1);
e_E = zeros(4,1);

for N = 1:4
    E_LCAO = zeros(size(R));
    E = zeros(size(R));
    for i = 1:length(R)
        S = exp(-dosRa0(i))*(1 + dosRa0(i) + dosRa0(i)^2/3);
        H11 = E1s + exp(-2*dosRa0(i))*(1/dosRa0(i) + 1)*alpha*hc/a0;
        H12 = 2*E1s*exp(-dosRa0(i))*(1/2 - 1/dosRa0(i) + 7/6*dosRa0(i) + (dosRa0(i)^2)/6);
        E_LCAO(i) = (H11 + H12)/(1 + S); 
        
        C = 2*me/(hc)^2*R(i)^2*(alpha*hc/2/R(i) - E_LCAO(i));
        B = 4*me/hc*alpha*R(i);
        if C <= 0
           C = 0.1;
        end 
        tol = 1e-12;
        max_iter = 100;
        h = 1e-8;
        for j = 1:max_iter
            frac = frac_cont(C,B,N);
            frac_h = frac_cont(C+h,B,N);
            diff = (frac_h - frac) / h;
            if diff == 0
                break;
            end
            C_nova = C - frac /diff;
            if abs(C_nova - C) < tol && abs(frac) < tol
                C = C_nova;
                break;
            end
            C = C_nova;
        end
        E(i) = alpha*hc/2/R(i) - hc^2*C/2/me/R(i)^2;
    end
    
    E_tots_ordres(N, :) = E;
    
    [E_min, idx_min] = min(E);
    R_min = dosRa0(idx_min);
    plot(dosRa0, E, 'Color', colors(N, :), ...
        'DisplayName', sprintf('Exacta (Matriu %dx%d). Mínim: (%.4fa_0, %.7f eV)', N,N, R_min, E_min));
    plot(R_min, E_min, 'o', 'MarkerEdgeColor', colors(N, :), ...
        'MarkerFaceColor', 'none', 'MarkerSize', 6, 'HandleVisibility', 'off');
    R_tots_ordres(N) = R_min;
    e_R(N) = max(dosRa0(idx_min-1)-R_min, dosRa0(idx_min+1)-R_min);
    e_E(N) = max(E(idx_min-1)-E_min, E(idx_min+1)-E_min);
end

[E_min_LCAO, idx_LCAO] = min(E_LCAO);
R_min_LCAO = dosRa0(idx_LCAO);
plot(dosRa0, E_LCAO, '--k', 'DisplayName', sprintf('LCAO. Mínim: (%.4fa_0, %.7f eV)', R_min_LCAO, E_min_LCAO));
plot([0 10], [E1s E1s], 'k-.', 'DisplayName', 'Energia E_{1s}');
plot(R_min_LCAO, E_min_LCAO, 'dk', 'MarkerFaceColor', 'none', 'MarkerSize', 6, 'HandleVisibility', 'off');
legend('Location', 'northeast');
ylim([-17, -12.5]); xlim([0, 10]);
xlabel('Distància internuclear (2R/a_0)'); ylabel('Energia (eV)');
title('Convergència del mètode exacte');
hold off;

%% 2. COMPARACIÓ AMB L'OSCIL·LADOR HARMÒNIC
[E_min, idx_min] = min(E);
dosRa0_min = dosRa0(idx_min);
h_deriv = dosRa0(2) - dosRa0(1);  
k = (E(idx_min+1) - 2*E(idx_min) + E(idx_min-1)) / h_deriv^2;
k_mes = (E(idx_min+2) - 2*E(idx_min+1) + E(idx_min)) / h_deriv^2;
k_menos = (E(idx_min) - 2*E(idx_min-1) + E(idx_min-2)) / h_deriv^2;
e_k = max(abs(k-k_mes),abs(k-k_menos));
oscil = E_min + 0.5*k*(dosRa0 - dosRa0_min).^2;

figure;
hold on;
plot(dosRa0, oscil, 'r--', 'DisplayName', sprintf('Oscil·lador harmònic, k = %.4f eV/a_0^2', k));
plot(dosRa0, E, 'b', 'DisplayName', 'Exacta (Matriu 4\times4)');
plot(dosRa0, E_LCAO, 'k--', 'DisplayName', 'LCAO');
plot([0 10], [E1s E1s], 'k-.', 'DisplayName', 'Energia E_{1s}');
plot(dosRa0_min, E_min, 'bo', 'MarkerFaceColor', 'none', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(R_min_LCAO, E_min_LCAO, 'dk', 'MarkerFaceColor', 'none', 'MarkerSize', 6, 'HandleVisibility', 'off');
xlabel('Distància internuclear (2R/a_0)'); ylabel('Energia (eV)');
title('Comparació amb l''aproximació d''oscil·lador harmònic');
legend('Location', 'northeast');
ylim([-17, -12.5]); xlim([0, 10]);
grid off; hold off;

%% 3. NUMEROV
mp = 938.27208816e6; % eV (29)
mu = mp/2; % eV
K = hc^2 / (2*mu); % eV·fm²

% Ajustem el límit superior a 6.0 per a permetre que les cues dels estats 
% excitats (v=1 i v=2) decandisquen correctament a zero abans de la paret de la caixa.
idx_retallada = find(dosRa0 >= 0.5 & dosRa0 <= 4.5); 
dosRa0_num = dosRa0(idx_retallada);
R_num = R(idx_retallada);
V_pot = E(idx_retallada)'; 
x = (R_num * 2)';  % 2R en fm per a Numerov
h_num = x(2) - x(1);                           
N_size = length(x); 

B = diag(10*ones(1, N_size)/12) + diag(ones(1, N_size-1)/12, 1) + diag(ones(1, N_size-1)/12, -1);
A = (diag(-2*ones(1, N_size)) + diag(ones(1, N_size-1), 1) + diag(ones(1, N_size-1), -1)) / h_num^2;
H = -K * (B \ A) + diag(V_pot);

[vectors, valors] = eig(H); 
En = real(diag(valors));
psi = vectors;
[En_ordenats, ordre] = sort(En);
psi_ordenats = psi(:, ordre);

% v=0
E_v0 = En_ordenats(1);
psi_v0 = psi_ordenats(:, 1);
psi_v0_norm = psi_v0 / sqrt(trapz(x, psi_v0.^2));

% v=1
E_v1 = En_ordenats(2);
psi_v1 = psi_ordenats(:, 2);
psi_v1_norm = psi_v1 / sqrt(trapz(x, psi_v1.^2));

% v=2
E_v2 = En_ordenats(3);
psi_v2 = psi_ordenats(:, 3);
psi_v2_norm = psi_v2 / sqrt(trapz(x, psi_v2.^2));


fprintf('\n========== ESTATS VIBRACIONALS (DOMINI RETALLAT) ==========\n');
fprintf('Mida efectiva de la matriu de Numerov: %d x %d\n', N_size, N_size);
fprintf('Energia v=0: %.7f eV\n', E_v0);
fprintf('Energia v=1: %.7f eV\n', E_v1);
fprintf('Energia v=2: %.7f eV\n', E_v2);

% ===== GRÀFICA PRINCIPAL AMB ELS TRES ESTATS =====
figure; hold on;
punts_tall0 = find(E <= E_v0);
punts_tall1 = find(E <= E_v1);
punts_tall2 = find(E <= E_v2);
c_v0 = [0.70 0.15 0.15]; c_v1 = [0.00 0.35 0.60]; c_v2 = [0.75 0.60 0.10]; 

fill([dosRa0(punts_tall2), fliplr(dosRa0(punts_tall2))], [E(punts_tall2), fliplr(E_v2 * ones(size(punts_tall2)))], c_v2, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([dosRa0(punts_tall1), fliplr(dosRa0(punts_tall1))], [E(punts_tall1), fliplr(E_v1 * ones(size(punts_tall1)))], c_v1, 'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([dosRa0(punts_tall0), fliplr(dosRa0(punts_tall0))], [E(punts_tall0), fliplr(E_v0 * ones(size(punts_tall0)))], c_v0, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

h_exc = plot(dosRa0, E, 'b-', 'LineWidth', 1.0);
h_lcao = plot(dosRa0, E_LCAO, 'k--');
h_h = plot([0 10], [E1s E1s], 'k-.');
h_min_exc = plot(R_min, E_min, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'none');
h_min_lcao = plot(R_min_LCAO, E_min_LCAO, 'kd', 'MarkerSize', 6);

h_v0 = line([dosRa0(punts_tall0(1)), dosRa0(punts_tall0(end))], [E_v0, E_v0], 'Color', c_v0, 'LineWidth', 1.0);
h_v1 = line([dosRa0(punts_tall1(1)), dosRa0(punts_tall1(end))], [E_v1, E_v1], 'Color', c_v1, 'LineWidth', 1.0);
h_v2 = line([dosRa0(punts_tall2(1)), dosRa0(punts_tall2(end))], [E_v2, E_v2], 'Color', c_v2, 'LineWidth', 1.0);

xlabel('Distància internuclear 2R/a_{0}'); ylabel('Energia (eV)');
title('Potencial electrònic i estats lligats de H_{2}^{+}');
legend([h_exc, h_lcao, h_h, h_min_exc, h_min_lcao, h_v0, h_v1, h_v2], ...
    {'Potencial exacte', 'Aproximació LCAO', 'Energia E_{1s}', ...
     sprintf('Mínim exacte: (%.4fa_0, %.7f eV)', R_min, E_min), ...
     sprintf('Mínim LCAO: (%.4fa_0, %.7f eV)', R_min_LCAO, E_min_LCAO), ...
     sprintf('v=0: %.7f eV', E_v0), sprintf('v=1: %.7f eV', E_v1), sprintf('v=2: %.7f eV', E_v2)}, 'Location', 'northeast');
xlim([0.5, 7]); ylim([-17, -12]); hold off;

% ===== GRÀFICA COMPARATIVA SUPERPOSADA (CORREGIDA AMB ELS NOUS LÍMITS) =====
figure; hold on;
xx_plot = x / a0;  
plot(xx_plot, psi_v0_norm, '-', 'Color', c_v0, 'LineWidth', 1.0, 'DisplayName', 'Fonamental (v=0)');
plot(xx_plot, psi_v1_norm, '-', 'Color', c_v1, 'LineWidth', 1.0, 'DisplayName', 'Excitat (v=1)');
plot(xx_plot, psi_v2_norm, '-', 'Color', c_v2, 'LineWidth', 1.0, 'DisplayName', 'Excitat (v=2)');
xlabel('Distància internuclear (2R/a_0)'); ylabel('\Psi');
title('Comparació de les funcions d''ona');
legend('Location', 'best'); 

% Sincronitzem el xlim estrictament amb els extrems reals de la retallada de dades
xlim([0.5, 4.5]); 
hold off;


%% ========================================================================
%  ESQUEMA DE NIVELLS D'ENERGIA ROTACIONALS (TAMBÉ AMB RETALLADA)
%  ========================================================================
v_valors = [0, 1, 2]; L_valors = [0, 1, 2, 3];
Energies_vL = zeros(length(v_valors), length(L_valors));
InvB_A = B \ A;

for l_idx = 1:length(L_valors)
    L_act = L_valors(l_idx);
    V_centrifug_act = K * L_act * (L_act+1) ./ (x.^2);
    V_eff_act = V_pot + V_centrifug_act;

    H_act = -K * InvB_A + diag(V_eff_act);
    [~, valors_act] = eig(H_act);
    En_act = sort(real(diag(valors_act)));

    for v_idx = 1:length(v_valors)
        Energies_vL(v_idx, l_idx) = En_act(v_valors(v_idx) + 1);
    end
end

fprintf('\n================ ESQUEMA D''ENERGIES (eV) ================\n');
fprintf('%-6s | %-12s | %-12s | %-12s | %-12s\n', 'v \ L', 'L = 0', 'L = 1', 'L = 2', 'L = 3');
fprintf('---------------------------------------------------------\n');
for v_idx = 1:length(v_valors)
    fprintf('v = %d   | %-12.7f | %-12.7f | %-12.7f | %-12.7f\n', ...
        v_valors(v_idx), Energies_vL(v_idx,1), Energies_vL(v_idx,2), Energies_vL(v_idx,3), Energies_vL(v_idx,4));
end

figure('Name', 'Esquema de nivells d''energia', 'NumberTitle', 'off'); hold on;
colors_L = [0 0 0; 0 0.2 0.6; 0.6 0 0; 0 0.5 0];     
x_esq = 0; x_dret = 1;  
for v_idx = 1:length(v_valors)
    v_act = v_valors(v_idx);
    for l_idx = 1:length(L_valors)
        E_nivell = Energies_vL(v_idx, l_idx);
        plot([x_esq, x_dret], [E_nivell, E_nivell], 'Color', colors_L(l_idx, :), 'LineWidth', 1.0);
    end
    E_min_bloc = Energies_vL(v_idx, 1);
    text(x_esq - 0.05, E_min_bloc, sprintf('v = %d', v_act), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
ylabel('Energia (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('Diagrama energètic de la molècula H_2^+');
xlim([-0.5, 1.5]); ylim([min(Energies_vL(:)) - 0.2, max(Energies_vL(:)) + 0.2]);
h_leg = zeros(1, 4);
for l_idx = 1:4, h_leg(l_idx) = plot(nan, nan, 'Color', colors_L(l_idx, :), 'LineWidth', 1.5); end
legend(h_leg, {'L = 0', 'L = 1', 'L = 2', 'L = 3'}, 'Location', 'northwest');
ax = gca; ax.XAxis.Visible = 'off'; hold off;


% CÀLCUL DE L'ERROR DE NUMEROV ANALÍTIC
h_adimensional = h_num / a0; 

% Ara la potència quarta actua com el percentatge de truncament pur exacte
error_v0 = (h_adimensional^4 / 12) * abs(E_v0);
error_v1 = (h_adimensional^4 / 12) * abs(E_v1);
error_v2 = (h_adimensional^4 / 12) * abs(E_v2);

fprintf('\n========== INCERTESES NUMÈRIQUES REALS CORREGIDES ==========\n');
fprintf('Error numèric real per a v=0: %.2e eV\n', error_v0);
fprintf('Error numèric real per a v=1: %.2e eV\n', error_v1);
fprintf('Error numèric real per a v=2: %.2e eV\n', error_v2);