%%%% Model 1 - Cobb Douglas, Exogeneous Human Capital %%%%
% L_a -> x(1), W_a -> x(2), W_m -> x(3)
%%%%%%

function F = Model1_Function(x, Z_at, Z_mt, alpha, Hbar_t, SD, P_t, T)
    F(1) = x(2) - (1 - alpha)*P_t * Z_at*(T^alpha) * x(1)^(-alpha); % wage a
    F(2) = x(3) - Z_mt; % wage m
    F(3) = x(1) - cdf('Normal',(x(2)/x(3)) - Hbar_t, 0 ,SD); % Labor supply
end