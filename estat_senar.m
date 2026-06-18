clear; clc; close all;

%% CONSTANTS I PARÀMETRES
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

dosRa0 = linspace(0.3, 10, 10000);
R_fm = dosRa0/2 * a0;  % R en fm

%% CÀLCUL LCAO
[E_LCAO_enll, S, H11, H12] = calcul_LCAO(dosRa0, alpha, hc, a0, E1s);
E_LCAO_anti = (H11 - H12) ./ (1 - S);      % estat antienllaçant (a=-1)

%% CÀLCUL SEPARACIÓ DE VARIABLES PER A N=1,2,3,4 (l senars)
colors = lines(4);
colors(3, :) = [0.10 0.55 0.25];

E_separacio_cells = cell(4, 1);
C_opt_cells = cell(4, 1);
info_cells = cell(4, 1);

tol = 1e-12;
max_iter = 100;
h_deriv = 1e-8;

for N = 1:4
    E_separacio = zeros(size(dosRa0));
    C_opt = zeros(size(dosRa0));
    
    for i = 1:length(dosRa0)
        % Estimació inicial de C a partir de l'energia LCAO antienllaçant
        C = 2 * me / hc^2 * R_fm(i)^2 * (alpha * hc / (2 * R_fm(i)) - E_LCAO_anti(i));
        B = 4 * me / hc * alpha * R_fm(i);
        
        % Si C és negativa o zero, posem un valor inicial petit positiu
        if C <= 0
            C = 0.1;
        end
        
        % Newton-Raphson per trobar C exacte (amb l senars)
        for j = 1:max_iter
            frac = frac_cont_imparell(C, B, N);
            frac_h = frac_cont_imparell(C + h_deriv, B, N);
            diff = (frac_h - frac) / h_deriv;
            
            if diff == 0
                break;
            end
            
            C_nova = C - frac / diff;
            
            if abs(C_nova - C) < tol && abs(frac) < tol
                C = C_nova;
                break;
            end
            
            C = C_nova;
        end
        
        C_opt(i) = C;
        E_separacio(i) = alpha * hc / (2 * R_fm(i)) - hc^2 * C / (2 * me * R_fm(i)^2);
    end
    
    % Guardar resultats
    E_separacio_cells{N} = E_separacio;
    C_opt_cells{N} = C_opt;
end

%% GRÀFICS DE CONVERGÈNCIA
figure('Name', 'Corbes d''energia de la funció electrònica senar', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);
hold on;

% Dibuixar per a N=1,2,3,4
for N_plot = 1:4
    E_plot = E_separacio_cells{N_plot};
    info_plot = info_cells{N_plot};
    
    plot(dosRa0, E_plot, 'Color', colors(N_plot, :), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('N=%d', N_plot));
end

% LCAO antienllaçant i E1s
plot(dosRa0, E_LCAO_anti, 'k--', 'LineWidth', 1.2, ...
     'DisplayName', sprintf('LCAO antienllaçant'));
plot([0 10], [E1s E1s], 'k-.', 'DisplayName', 'Energia E_{1s}');

xlabel('Distància internuclear (2R/a_0)', 'FontSize', 12);
ylabel('Energia (eV)', 'FontSize', 12);
title('Convergència del mètode de separació de variables (l senars)', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 10); 
ylim([-14, 0]); xlim([0, 10]);
grid on;
hold off;


%% FUNCIÓ AUXILIAR

function frac = frac_cont_imparell(C, B, N)
    % Fracció contínua per a l senars
    if C <= 0
        C = 1e-10;
    end
    
    s = B/(2*sqrt(C)) - 1;
    
    % Càlcul de la constant de separació A per a l senars
    switch N
        case 1
            % Només l = 1
            A = 2 - (2/5)*C;
            
        case 2
            % l = 1, 3
            M = zeros(2);
            % l = 1
            M(1,1) = -2 + (2/5)*C;
            % l = 3
            M(2,2) = -12 + (36/55)*C;
            % Acoblament
            M(1,2) = (4/15)*C;
            M(2,1) = (12/35)*C;
            valors_propis = eig(-M);
            A = min(valors_propis);
            
        case 3
            % l = 1, 3, 5
            M = zeros(3);
            ll = [1, 3, 5];
            for idx = 1:3
                l = ll(idx);
                M(idx,idx) = -l*(l+1) + C*((l+1)^2/((2*l+1)*(2*l+3)) + l^2/(4*l^2-1));
                if idx > 1
                    l_ant = ll(idx-1);
                    M(idx, idx-1) = C*(l_ant*(l_ant+1))/((2*l_ant+1)*(2*l_ant+3));
                end
                if idx < 3
                    l_seg = ll(idx);
                    M(idx, idx+1) = C*((l_seg+2)*(l_seg+1))/((2*l_seg+5)*(2*l_seg+3));
                end
            end
            valors_propis = eig(-M);
            A = min(valors_propis);
            
        case 4
            % l = 1, 3, 5, 7
            M = zeros(4);
            ll = [1, 3, 5, 7];
            for idx = 1:4
                l = ll(idx);
                M(idx,idx) = -l*(l+1) + C*((l+1)^2/((2*l+1)*(2*l+3)) + l^2/(4*l^2-1));
                if idx > 1
                    l_ant = ll(idx-1);
                    M(idx, idx-1) = C*(l_ant*(l_ant+1))/((2*l_ant+1)*(2*l_ant+3));
                end
                if idx < 4
                    l_seg = ll(idx);
                    M(idx, idx+1) = C*((l_seg+2)*(l_seg+1))/((2*l_seg+5)*(2*l_seg+3));
                end
            end
            valors_propis = eig(-M);
            A = min(valors_propis);
    end
    
    % Fracció contínua
    b = @(n) -(2*n.^2 + (4*sqrt(C) - 2*s)*n + C + A - 2*sqrt(C)*s - s) / (n+1)^2;
    a = @(n) -(n - 1 - s).^2 / (n+1)^2;
    
    tiny = 1e-30;
    f = b(0);
    if f == 0
        f = tiny;
    end
    c = f;
    d = 0;
    
    for n = 1:200
        dj = b(n) + a(n)*d;
        if dj == 0
            dj = tiny;
        end
        cj = b(n) + a(n)/c;
        if cj == 0
            cj = tiny;
        end
        d = 1/dj;
        c = cj;
        dif = c*d;
        f = f*dif;
        if abs(dif - 1) < 1e-15
            break;
        end
    end
    frac = f;
end