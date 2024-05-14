%%%% Final Model %%%%
% L_a -> x(1), W_a -> x(2), W_m -> x(3)
%%%%%%

function F = Final_Model_Function(x, Z_lt, Z_tt, Z_mt, Mu, Hbar_t, SD, P_t, T_t)
    F(1) = x(2) - P_t*Z_lt*(1 + ((Z_tt*T_t)/(Z_lt*x(1)))^((Mu-1)/Mu))...
        ^(1/(Mu - 1)); % wage a
    F(2) = x(3) - Z_mt; % wage m
    F(3) = x(1) - cdf('Normal',(x(2)/x(3)) - Hbar_t, 0 ,SD); % Labor supply
end