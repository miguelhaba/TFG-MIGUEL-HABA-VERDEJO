%% ========================================================================
%  GRAFICACIÓ DE LA FUNCIÓ D'ONA EXACTA DE L'ELECTRÓ (BATES / JAFFÉ / HYLLERAAS)
%  ========================================================================
clear; clc; close all;
% --- 1. Constants Físiques i Paràmetres de l'estat d'equilibri (H2+) ---
alpha = 7.2973525693e-3; 
hc = 197.3269804e6; % eV·fm
me = 0.51099895e6;  % eV
a0 = hc / (me * alpha); % fm
% Prenem els valors clàssics d'equilibri obtinguts del mètode exacte
% per a evitar dependències de funcions externes:
dosRa0_eq = 1.9977;          % Distància internuclear d'equilibri (2R/a0)
R_eq = dosRa0_eq / 2;      % R en unitats d'a0
E_total = -16.3985231;        % Energia total calculada en eV
E1s = -alpha * hc / 2 / a0; % Energia de l'àtom d'H (~ -13.6057 eV)
% Conversió de l'energia electrònica a unitats atòmiques (Rydbergs/Hartrees equivalents)
E_elec = E_total - (alpha * hc / 2 / R_eq / a0); % Restem repulsió nuclear
C_param = (R_eq * a0)^2 * (-2 * me * E_elec / hc^2); % El paràmetre C adimensional
B_param = 4 * me / hc * alpha * (R_eq * a0);         % El paràmetre B adimensional
% --- 2. Càlcul dels coeficients de l'equació Angular g(eta) ---
% Resolem el problema de valors propis per a la matriu M (ordre 6x6)
N_angular = 4;
M = zeros(N_angular);
l_vals = 0:2:(2*N_angular-2); % Per a l'estat fonamental, només l parells per simetria g
for idx = 1:N_angular
    l = l_vals(idx);
    % Terme diagonal (sense A0)
    M(idx, idx) = -l*(l+1) + C_param * ((l+1)^2 / ((2*l+1)*(2*l+3)) + l^2 / (4*l^2-1));
    % Termes de la sub i superdiagonal
    if idx > 1
        l_prev = l_vals(idx-1);
        M(idx, idx-1) = C_param * (l*(l-1)) / ((2*l-3)*(2*l-1));
    end
    if idx < N_angular
        l_next = l_vals(idx+1);
        M(idx, idx+1) = C_param * ((l+2)*(l+1)) / ((2*l+5)*(2*l+3));
    end
end
% Els valors propis de -M corresponen a les constants de separació A0
[V_ang, D_ang] = eig(-M);
[A0_val, idx_A0] = min(diag(D_ang)); % El valor propi més baix és l'estat fonamental
c_coefs = V_ang(:, idx_A0);         % Coeficients c_l angulars
% --- 3. Càlcul dels coeficients de l'equació Radial f(xi) ---
N_radial = 8;
sigma = (B_param / (2 * sqrt(C_param))) - 1;
A = zeros(N_radial);
for n = 0:N_radial-1
    idx = n + 1;
    % Coeficient alfa_n multiplicat per (n+1)^2 per a fer la matriu del sistema lineal
    A(idx, idx) = -(2*n^2 + (4*sqrt(C_param) - 2*sigma)*n + C_param + A0_val - 2*sigma*sqrt(C_param) - sigma);
    if idx < N_radial
        A(idx, idx+1) = (n+1)^2;
    end
    if idx > 1
        A(idx, idx-1) = (n - sigma - 1)^2;
    end
end
% Busquem el nucli de la matriu A (sistema homogeni A*d = 0)
[V_rad, ~] = eig(A);
d_coefs = V_rad(:, 1); % Coeficients d_n radials (sèrie de Jaffé)
% --- 4. Construcció de la Malla Espacial en el Pla XZ ---
x_mesh = linspace(-5, 5, 500);
z_mesh = linspace(-5, 5, 500);
[X, Z] = meshgrid(x_mesh, z_mesh);
% Posicions reals dels nuclis (0, +R) i (0, -R)
nucli1_x = 0; nucli1_z =  R_eq;
nucli2_x = 0; nucli2_z = -R_eq;
% Distàncies atòmiques r1 i r2 des de cada punt de la malla
r1 = sqrt((X - nucli1_x).^2 + (Z - nucli1_z).^2);
r2 = sqrt((X - nucli2_x).^2 + (Z - nucli2_z).^2);
% Transformació a coordenades el·líptiques/esferoïdals prolates (xi, eta)
Xi = (r1 + r2) / (2 * R_eq);
Eta = (r1 - r2) / (2 * R_eq);
% Estabilització de contorns per a evitar divergències numèriques flotants
Xi(Xi < 1) = 1;
Eta(Eta > 1) = 1; Eta(Eta < -1) = -1;
% Variable de Jaffé u
U_var = (Xi - 1) ./ (Xi + 1);
% --- 5. Avaluació Numèrica de les Funcions d'Ona Exactes ---
F_xi = zeros(size(X));
G_eta = zeros(size(X));
% Avaluació de la part radial f(xi) combinant la sèrie de potències
omega_u = zeros(size(X));
for n = 0:N_radial-1
    omega_u = omega_u + d_coefs(n+1) * (U_var.^n);
end
F_xi = exp(-sqrt(C_param) * Xi) .* ((1 - U_var).^(-sigma)) .* omega_u;
% Avaluació de la part angular g(eta) combinant els Polinomis de Legendre
for idx = 1:N_angular
    l = l_vals(idx);
    % Generació del polinomi de Legendre d'ordre l de manera matricial
    P_l = legendre(l, Eta);
    % legendre() retorna una matriu on la primera fila és m=0
    P_l0 = squeeze(P_l(1, :, :)); 
    if l == 0, P_l0 = ones(size(Eta)); end % Cas base estable
    
    G_eta = G_eta + c_coefs(idx) * P_l0;
end
% Funció d'ona exacta total (Producte de les dues parts)
Psi_exacta = F_xi .* G_eta;
% Normalització visual del pic d'intensitat sobre el valor formal de l'orbital 1s
Psi_exacta = Psi_exacta / max(Psi_exacta(:)) * (1/sqrt(pi));
% --- 6. Generació dels Gràfics de Contorns 2D ---
figure('Position', [150, 150, 1100, 500], 'Name', 'Funció d''ona exacta H2+', 'NumberTitle', 'off');
contourf(X, Z, Psi_exacta, 20);
colorbar;
hold on;
xlabel('x (a_0)');
ylabel('z (a_0)');
title(sprintf('Funció d''ona per a 2R = %.4fa_0', dosRa0_eq));
axis equal;
colormap(jet);