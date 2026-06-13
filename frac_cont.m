%% Algorisme de Lentz modificat

function frac = frac_cont(C,B,N)
if C <= 0
   C = 1e-10;
end 

s = B/(2*sqrt(C)) - 1;
switch N
    case 1  
        A = -C/3;
    case 2
        A = -(3/7*C-3) - 3*sqrt(1 - 4/63*C + 8/735*C^2);
    case 3
        M = [-C/3, -(2/15)*C, 0; -(2/3)*C, 6 - (11/21)*C, -(4/21)*C; 0, -(12/35)*C, 20 - (39/77)*C];
        valors_propis = eig(M);
        A = min(valors_propis);
    case 4
        % Construir la matriu amb la recurrència
        M = zeros(4,4);
        ll = [0, 2, 4, 6];
        for idx = 1:4
            l = ll(idx);
            M(idx,idx) = -l*(l+1) + C*((l+1)^2/((2*l+1)*(2*l+3)) + l^2/(4*l^2-1));
            if idx > 1
                l_ant = ll(idx);
                M(idx, idx-1) = C*(l_ant*(l_ant-1))/((2*l_ant-3)*(2*l_ant-1));
            end
            if idx < 4
                interlayer = l;
                M(idx, idx+1) = C*((interlayer+2)*(interlayer+1))/((2*interlayer+5)*(2*interlayer+3));
            end
        end
        valors_propis = eig(-M); 
        A = min(valors_propis);
end

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
    if n == 200
        fprintf('S''ha arribat al màxim d''iteracions.')
    end
end
frac = f;
end