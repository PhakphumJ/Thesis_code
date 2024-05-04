%%%% Model 1 - Cobb Douglas, Exogeneous Human Capital %%%%
% L_a -> x(1), Y_a -> x(2), Y_m -> x(3), W_a -> x(4), W_m -> x(5)
%%%%%%

function F = Model1_Function(x, Z_at, Z_mt, alpha, Hbar_t, SD, P_t, T)
    F(1) = x(2) - Z_at*(T^alpha) * x(1)^(1-alpha); % Production function a
    F(2) = x(3) - Z_mt*(1 - x(1)); % Production function m
    F(3) = x(4) - (1 - alpha)*P_t * Z_at*(T^alpha) * x(1)^(-alpha); % wage a
    F(4) = x(5) - Z_mt; % wage m
    F(5) = x(1) - cdf('Normal',(x(4)/x(5)) - Hbar_t, 0 ,SD); % Labor supply
end