function resultats = numerov(dosRa0, R_fm, E_potencial, a0, K)

% Retallar el domini
idx_retallada = find(dosRa0 >= 0.5 & dosRa0 <= 4.5);
dosRa0_num = dosRa0(idx_retallada);
R_num = R_fm(idx_retallada);
V_pot = E_potencial(idx_retallada)';

x = (R_num * 2)';  % 2R en fm
h_num = x(2) - x(1);
N_size = length(x);

% Matrius de Numerov
B = diag(10*ones(1, N_size)/12) + diag(ones(1, N_size-1)/12, 1) + diag(ones(1, N_size-1)/12, -1);
A_mat = (diag(-2*ones(1, N_size)) + diag(ones(1, N_size-1), 1) + diag(ones(1, N_size-1), -1)) / h_num^2;

% Hamiltonià
H = -K * (B \ A_mat) + diag(V_pot);

% Diagonalització
[vectors, valors] = eig(H);
En = real(diag(valors));
psi = vectors;
[En_ordenats, ordre] = sort(En);
psi_ordenats = psi(:, ordre);

% Extraure estats v=0,1,2
resultats.v = cell(3, 1);
for iv = 0:2
    resultats.v{iv+1}.E = En_ordenats(iv+1);
    psi_norm = psi_ordenats(:, iv+1);
    resultats.v{iv+1}.psi = psi_norm / sqrt(trapz(x, psi_norm.^2));
end

% Nivells rotacionals
L_valors = [0, 1, 2, 3];
resultats.Energies_vL = zeros(3, length(L_valors));
InvB_A = B \ A_mat;

for l_idx = 1:length(L_valors)
    L_act = L_valors(l_idx);
    V_centrifug = K * L_act * (L_act+1) ./ (x.^2);
    V_eff = V_pot + V_centrifug;

    H_act = -K * InvB_A + diag(V_eff);
    [~, valors_act] = eig(H_act);
    En_act = sort(real(diag(valors_act)));

    for v_idx = 1:3
        resultats.Energies_vL(v_idx, l_idx) = En_act(v_idx);
    end
end

resultats.x = x;
resultats.dosRa0_num = dosRa0_num;
resultats.V_pot = V_pot;
resultats.L_valors = L_valors;
resultats.h_num = h_num;
resultats.N_size = N_size;
end