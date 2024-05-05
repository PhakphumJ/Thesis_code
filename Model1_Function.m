%%%% Model 1 - Cobb Douglas, Exogeneous Human Capital %%%%
% L_a -> x(1), W_a -> x(2), W_m -> x(3)
%%%%%%

function F = Model1_Function(x, Z_at, Z_mt, alpha, Hbar_t, SD, P_t, T)
    F(1) = x(2) - (1 - alpha)*P_t * Z_at*(T^alpha) * x(1)^(-alpha); % wage a
    F(2) = x(3) - Z_mt; % wage m
    F(3) = x(1) - cdf('Normal',(x(2)/x(3)) - Hbar_t, 0 ,SD); % Labor supply
    F(4) = x(4) - log(x(1)); % Log labor supply
    F(5) = x(5) - (1/alpha)*log(1-alpha) - log(T) - (1/alpha)*(log(Z_at/Z_mt) + log(P_t))...
        + (1/alpha)*log(x(2)/x(3)); % Log labor demand
    F(6) = x(5) - x(4); % log labor demand = log labor supply

    % I am not sure if I need F(4), F(5), and F(6). 
    % Adding it does not change the results when 
    % I tried using it in the file "Model1_OnePeriod_Try.m"
end