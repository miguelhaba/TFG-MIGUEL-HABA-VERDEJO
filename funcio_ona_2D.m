function [Psi_exacta, X, Z, info_ona] = funcio_ona_2D(R_eq, C_param, B_param, N_angular, N_radial)
    
    if nargin < 6, N_angular = 4; end
    if nargin < 7, N_radial = 8; end
    
% Càlcul dels coeficients angulars
l_vals = 0:2:(2*N_angular-2);
M = zeros(N_angular);
    
for idx = 1:N_angular
    l = l_vals(idx);
    M(idx, idx) = -l*(l+1) + C_param * ((l+1)^2/((2*l+1)*(2*l+3)) + l^2/(4*l^2-1));
    if idx > 1
       M(idx, idx-1) = C_param * (l*(l-1))/((2*l-3)*(2*l-1));
    end
    if idx < N_angular
       M(idx, idx+1) = C_param * ((l+2)*(l+1))/((2*l+5)*(2*l+3));
    end
end
    
[V_ang, D_ang] = eig(-M);
[A0_val, idx_A0] = min(diag(D_ang));
c_coefs = V_ang(:, idx_A0);
    
% Càlcul dels coeficients radials
sigma = (B_param / (2 * sqrt(C_param))) - 1;
A_mat = zeros(N_radial);
    
for n = 0:N_radial-1
    idx = n + 1;
    A_mat(idx, idx) = -(2*n^2 + (4*sqrt(C_param) - 2*sigma)*n + C_param + A0_val - 2*sigma*sqrt(C_param) - sigma);
    if idx < N_radial
       A_mat(idx, idx+1) = (n+1)^2;
    end
    if idx > 1
       A_mat(idx, idx-1) = (n - sigma - 1)^2;
    end
end
    
[V_rad, ~] = eig(A_mat);
d_coefs = V_rad(:, 1);
    
% Malla 
x_mesh = linspace(-5, 5, 500);
z_mesh = linspace(-5, 5, 500);
[X, Z] = meshgrid(x_mesh, z_mesh);
    
% Posicions dels nuclis
nucli1_x = 0; nucli1_z =  R_eq;
nucli2_x = 0; nucli2_z = -R_eq;

% Calcular distàncies als nuclis
r1 = sqrt((X - nucli1_x).^2 + (Z - nucli1_z).^2);
r2 = sqrt((X - nucli2_x).^2 + (Z - nucli2_z).^2);

% Funció d'ona per separació de variables
Xi = (r1 + r2) / (2 * R_eq);
Eta = (r1 - r2) / (2 * R_eq);
  
Xi(Xi < 1) = 1;
Eta(Eta > 1) = 1; Eta(Eta < -1) = -1;
    
U_var = (Xi - 1) ./ (Xi + 1);
   
% Part radial
omega_u = zeros(size(X));
for n = 0:N_radial-1
    omega_u = omega_u + d_coefs(n+1) * (U_var.^n);
end
F_xi = exp(-sqrt(C_param) * Xi) .* ((1 - U_var).^(-sigma)) .* omega_u;
    
% Part angular
G_eta = zeros(size(X));
for idx = 1:N_angular
    l = l_vals(idx);
    P_l = legendre(l, Eta);
    P_l0 = squeeze(P_l(1, :, :));
    if l == 0, P_l0 = ones(size(Eta)); end
    G_eta = G_eta + c_coefs(idx) * P_l0;
end
    
Psi_exacta = F_xi .* G_eta;
Psi_exacta = abs(Psi_exacta / max(abs(Psi_exacta(:))) * (1/sqrt(pi)));
    
% Guardar informació
info_ona.c_coefs = c_coefs;
info_ona.d_coefs = d_coefs;
info_ona.A0_val = A0_val;
info_ona.l_vals = l_vals;
info_ona.sigma = sigma;
info_ona.N_radial = N_radial;
info_ona.N_angular = N_angular;
end