%% Mètode LCAO
clear; clc; close all;
% Constants
alpha = 7.2973525693e-3; %   ERROR (11)
c = 299792458e15; % fm/s  EXACTE  
h = 6.62607015e-34*1.602176634e-19/2/pi; % eV·s  EXACTE
hc = h*c; % eV·fm  EXACTE
me = 0.51099895000e6; %eV  ERROR (15)
a0 = hc/(me*alpha); % fm  
E1s = -alpha*hc/2/a0; % eV
dosRa0 = linspace(0.3,10,10000);


S = exp(-dosRa0).*(1 + dosRa0 + (dosRa0.^2)/3);
H11 = E1s + exp(-2*dosRa0).*(1./dosRa0 + 1)*alpha*hc/a0;
H12 = 2*E1s*exp(-dosRa0).*(1/2 - 1./dosRa0 + 7/6*dosRa0 + (dosRa0.^2)/6);

a = [1 -1];

figure
hold on
for i = 1:2
    E = (H11 + a(i)*H12)./(1 + a(i)*S);
    plot(dosRa0,E);
end
plot([0 10],[E1s E1s],'k--')
ylim([-16, 0]);
xlabel('Distància internuclear (2R / a_0)');
ylabel('Energia (eV)');
title('Mètode LCAO: Fites per a l''energia de l''estat estacionari');
legend('a = 1 (Enllaçant)', 'a = -1 (Antienllaçant)', 'Energia E_{1s}', 'Location', 'best');
hold off

% Gràfic només amb l'estat enllaçant
E = (H11 + H12)./(1 + S);
[E_min, idx_min] = min(E);
dosR_min = dosRa0(idx_min);
figure
hold on
plot(dosRa0,E,'k');
plot(dosRa0(idx_min),E_min,'ko');
plot([0 10],[E1s E1s],'k--')
ylim([-16, 0]);
xlabel('Distància internuclear (2R / a_0)');
ylabel('Energia (eV)');
title('Mètode LCAO: Fita per a l''energia de l''estat estacionari');
text_minim = sprintf('Mínim (2R = %.4f a_0, E = %.7f eV)', dosRa0(idx_min), E_min);
legend('Fita LCAO', text_minim, 'Energia E_{1s}', 'Location', 'best');
hold off

% Error en el mínim
e_R = max(dosRa0(idx_min-1)-dosR_min, dosRa0(idx_min+1)-dosR_min);
e_E = max(E(idx_min-1)-E_min, E(idx_min+1)-E_min);

%% Gràfic PSI
dosRa0 = dosRa0(idx_min);  % 2R en unitats de a_0
R = dosRa0/2;     % R en unitats de a_0

% Solapament per a 2R = 2.494 a_0
S_val = exp(-dosRa0)*(1 + dosRa0 + dosRa0^2/3);

% Crear malla de punts en el pla xz (tall transversal)
x = linspace(-5, 5, 500);
z = linspace(-5, 5, 500);
[X, Z] = meshgrid(x, z);

% Posicions dels nuclis: (0, 0, +R) i (0, 0, -R)
% En el tall xz, el nucli 1 està en (0, +R) i el nucli 2 en (0, -R)
nucli1_x = 0;
nucli1_z = R;
nucli2_x = 0;
nucli2_z = -R;

% Distàncies als nuclis en el pla xz
r1 = sqrt((X - nucli1_x).^2 + (Z - nucli1_z).^2);
r2 = sqrt((X - nucli2_x).^2 + (Z - nucli2_z).^2);

% Funcions d'ona atòmiques (1s de l'hidrogen) en unitats de a_0
chi1 = 1/sqrt(pi) * exp(-r1);
chi2 = 1/sqrt(pi) * exp(-r2);

% Estat enllaçant i antienllaçant NORMALITZATS amb 1/sqrt(2(1±S))
Psi_enllacant = (chi1 + chi2) / sqrt(2*(1 + S_val));
Psi_antienllacant = (chi1 - chi2) / sqrt(2*(1 - S_val));

%% Gràfics 2D en el pla xz
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Estat enllaçant
subplot(1, 2, 1);
contourf(X, Z, Psi_enllacant, 20);
colorbar;
hold on;
xlabel('x (a_0)');
ylabel('z (a_0)');
title(sprintf('Funció d''ona enllaçant per a 2R = %.4f a_0', dosRa0));
axis equal;
colormap(jet);

% Subplot 2: Estat antienllaçant
subplot(1, 2, 2);
contourf(X, Z, Psi_antienllacant, 20);
colorbar;
hold on;
xlabel('x (a_0)');
ylabel('z (a_0)');
title(sprintf('Funció d''ona antienllaçant per a 2R = %.4f a_0', dosRa0));
axis equal;
colormap(jet);

%% Tall 1D al llarg de l'eix internuclear (eix z, x=0)
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
xlabel('z (a_0)');
ylabel('\Psi');
title(sprintf('Tall de la funció d''ona al llarg de l''eix internuclear (2R = %.4f a_0)', dosRa0));
legend('Enllaçant', 'Antienllaçant', 'Àtom aïllat, \chi_1', 'Location', 'best');
grid off;