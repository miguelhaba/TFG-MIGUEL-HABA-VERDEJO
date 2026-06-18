function [E_exacta, C_opt, info] = calcul_separacio_variables(dosRa0, R_fm, E_LCAO, alpha, hc, me, N)
E_exacta = zeros(size(dosRa0));
C_opt = zeros(size(dosRa0));

tol = 1e-12;
max_iter = 100;
h = 1e-8;

for i = 1:length(dosRa0)
    % Valor inicial de C a partir de LCAO
    C = 2 * me / hc^2 * R_fm(i)^2 * (alpha * hc / (2 * R_fm(i)) - E_LCAO(i));
    B = 4 * me / hc * alpha * R_fm(i);

    if C <= 0
        C = 0.1;
    end

    % Newton-Raphson per trobar C exacte
    for j = 1:max_iter
        frac = frac_cont(C, B, N);
        frac_h = frac_cont(C + h, B, N);
        diff = (frac_h - frac) / h;

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
    E_exacta(i) = alpha * hc / (2 * R_fm(i)) - hc^2 * C / (2 * me * R_fm(i)^2);
end

% Informació del mínim
[E_min, idx_min] = min(E_exacta);
info.E_min = E_min;
info.idx_min = idx_min;
info.dosRa0_min = dosRa0(idx_min);
info.R_min = info.dosRa0_min / 2; % en unitats a0
info.C_min = C_opt(idx_min);
info.B_min = 4 * me / hc * alpha * R_fm(idx_min);
end