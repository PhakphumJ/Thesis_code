%%%% Model 2 - CES, Exogeneous Human Capital %%%%
% L_a -> x(1), W_a -> x(2), W_m -> x(3)
%%%%%%

function F = Model1_Function(x, Z_lt, Z_tt, Z_mt, Mu, Gamma, Hbar_t, SD, P_t, T)
    F(1) = x(2) - P_t*Z_lt*Gamma*(Gamma + (1 - Gamma)*(Z_tt*T/(Z_lt*x(1))...
        )^((Mu-1)/Mu))^(1/(Mu - 1)); % wage a
    F(2) = x(3) - Z_mt; % wage m
    F(3) = x(1) - cdf('Normal',(x(2)/x(3)) - Hbar_t, 0 ,SD); % Labor supply
end