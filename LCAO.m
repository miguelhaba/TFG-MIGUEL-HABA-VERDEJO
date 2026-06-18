clear; clc; close all;

% Constants
alpha = 7.2973525693e-3;
c = 299792458e15; % fm/s
h = 6.62607015e-34 * 1.602176634e-19 / (2*pi); % eV·s
hc = h * c; % eV·fm
me = 0.51099895000e6; % eV
a0 = hc / (me * alpha); % fm
E1s = -alpha * hc / (2 * a0); % eV

%% CÀLCUL
dosRa0 = linspace(0.3, 10, 10000);

% Calcular integrals LCAO
[E_LCAO, S, H11, H12] = calcul_LCAO(dosRa0, alpha, hc, a0, E1s);

% Estat antienllaçant
E_anti = (H11 - H12) ./ (1 - S);

% Trobar el mínim de l'estat enllaçant
[E_min, idx_min] = min(E_LCAO);
dosRa0_min = dosRa0(idx_min);
R_min = dosRa0_min / 2;

fprintf('========== RESULTATS LCAO ==========\n');
fprintf('Distància d''equilibri: 2R = %.4f a0\n', dosRa0_min);
fprintf('Energia mínima: E = %.7f eV\n', E_min);
fprintf('====================================\n\n');


%% GRÀFICS
% Gràfic estats enllaçant i antienllaçant
figure
hold on
plot(dosRa0, E_LCAO, 'b-', 'LineWidth', 1.5);
plot(dosRa0, E_anti, 'r-', 'LineWidth', 1.5);
plot([0 10], [E1s E1s], 'k--', 'LineWidth', 1);
ylim([-16, 0]);
xlabel('Distància internuclear (2R / a_0)', 'FontSize', 12);
ylabel('Energia (eV)', 'FontSize', 12);
title('Mètode LCAO: Fites per a l''energia de l''estat estacionari', 'FontSize', 14);
legend('a = 1 (Enllaçant)', 'a = -1 (Antienllaçant)', 'Energia E_{1s}', 'Location', 'best', 'FontSize', 12);
hold off

% Gràfic 2D
x = linspace(-5, 5, 500);
z = linspace(-5, 5, 500);
[X, Z] = meshgrid(x, z);

% Posicions dels nuclis
nucli1_x = 0; nucli1_z = R_min;
nucli2_x = 0; nucli2_z = -R_min;

% Distàncies als nuclis
r1 = sqrt((X - nucli1_x).^2 + (Z - nucli1_z).^2);
r2 = sqrt((X - nucli2_x).^2 + (Z - nucli2_z).^2);

% Funcions d'ona atòmiques
chi1 = 1/sqrt(pi) * exp(-r1);
chi2 = 1/sqrt(pi) * exp(-r2);

% Solapament a la distància d'equilibri
S_val = exp(-dosRa0_min) * (1 + dosRa0_min + dosRa0_min^2/3);

% Estats enllaçant i antienllaçant normalitzats
Psi_enllacant = (chi1 + chi2) / sqrt(2*(1 + S_val));
Psi_antienllacant = (chi1 - chi2) / sqrt(2*(1 - S_val));

figure('Position', [100, 100, 1200, 500]);

% Estat enllaçant
subplot(1, 2, 1);
contourf(X, Z, Psi_enllacant, 20);
colorbar;
hold on;
xlabel('x (a_0)', 'FontSize', 12);
ylabel('z (a_0)', 'FontSize', 12);
title(sprintf('Funció d''ona enllaçant per a 2R = %.4f a_0', dosRa0_min), 'FontSize', 14);
axis equal;
colormap(jet);

% Estat antienllaçant
subplot(1, 2, 2);
contourf(X, Z, Psi_antienllacant, 20);
colorbar;
hold on;
xlabel('x (a_0)', 'FontSize', 12);
ylabel('z (a_0)', 'FontSize', 12);
title(sprintf('Funció d''ona antienllaçant per a 2R = %.4f a_0', dosRa0_min), 'FontSize', 14);
axis equal;
colormap(jet);

% Gràfic 1D 
z_1d = linspace(-5, 5, 1000);
r1_1d = abs(z_1d - nucli1_z);
r2_1d = abs(z_1d - nucli2_z);

chi1_1d = 1/sqrt(pi) * exp(-r1_1d);
chi2_1d = 1/sqrt(pi) * exp(-r2_1d);

Psi_enll_1d = (chi1_1d + chi2_1d) / sqrt(2*(1 + S_val));
Psi_anti_1d = (chi1_1d - chi2_1d) / sqrt(2*(1 - S_val));

figure;
plot(z_1d, Psi_enll_1d, 'b-', 'LineWidth', 1.5);
hold on;
plot(z_1d, Psi_anti_1d, 'r-', 'LineWidth', 1.5);
plot(z_1d, chi1_1d, 'k:', 'LineWidth', 1);
xline(nucli1_z, '--k', 'Nucli 1');
xline(nucli2_z, '--k', 'Nucli 2');
xlabel('z (a_0)', 'FontSize', 12);
ylabel('\phi', 'FontSize', 12);
title(sprintf('Tall de la funció d''ona al llarg de l''eix internuclear (2R = %.4f a_0)', dosRa0_min), 'FontSize', 14);
legend('Enllaçant', 'Antienllaçant', 'Àtom aïllat, \chi_1', 'Location', 'best', 'FontSize', 12);
grid off;