cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
%% Declare parameters values
% Production parameters
Z_l_0 = 30;
Z_t_0 = 10;
Z_m_0 = 100;
Gamma = 0.5;
Mu = 0.5;

g_z_l = 0.04;
g_z_t = 0.03;
g_z_m = 0.02;

% Human capital parameters
H_bar_0 = 0.1;
g_H_bar = 0.1;
SD = 1;

% Price
P_0 = 2;
g_P = 0.01;


N = 30;
%% Create Exogeneous Variables (Z, p, H_bar, SD)
Z_l = [];
Z_l = [Z_l; Z_l_0];
for t = 1:N-1
   z_l_t = Z_l(t);
   Z_l_t_plus1 = z_l_t*(1 + g_z_l);
   Z_l = [Z_l; Z_l_t_plus1];
end

Z_t = [];
Z_t = [Z_t; Z_t_0];
for t = 1:N-1
   z_t_t = Z_t(t);
   Z_t_t_plus1 = z_t_t*(1 + g_z_t);
   Z_t = [Z_t; Z_t_t_plus1];
end

Z_m = [];
Z_m = [Z_m; Z_m_0];
for t = 1:N-1
   z_m_t = Z_m(t);
   Z_m_t_plus1 = z_m_t*(1 + g_z_m);
   Z_m = [Z_m; Z_m_t_plus1];
end


P = [];
P = [P; P_0];
for t = 1:N-1
    p_t = P(t);
    p_t_plus1 = p_t*(1 + g_P);
    P = [P; p_t_plus1];
end


H_bar = [];
H_bar = [H_bar; H_bar_0];
for t = 1:N-1
    H_bar_t = H_bar(t);
    H_bar_t_plus1 = H_bar_t*(1 + g_H_bar);
    H_bar = [H_bar; H_bar_t_plus1];
end


T = 10; %Land
%% Solve the model for each t
% Description: L_a -> x(1), W_a -> x(2), W_m -> x(3)

x0 = [0.5; 20; 2]; % Just a starting point for fsolve.

endo_mat = zeros(N,3); % To store results from all t.
Real_world_W_m_vec = zeros(N,1); %To store results from all t.
Y_a_t_mat = zeros(N,1); %To store results from all t.
Y_m_t_mat = zeros(N,1); %To store results from all t.
Output_perwk_vec = zeros(N,1); %To store results from all t.

options = optimoptions(@fsolve, 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, 'FunctionTolerance', 1.0e-05);

for t = 1:N
    Z_lt = Z_l(t);
    Z_tt = Z_t(t);
    Z_mt = Z_m(t);
    Hbar_t = H_bar(t);
    P_t = P(t);

    % Find the equilibrium
    x_star = fsolve(@(x)Model2_Function(x, Z_lt, Z_tt, Z_mt, Mu, Gamma, Hbar_t, SD, P_t, T), x0);
    
    % Calculate W_m * avg h in sector m
    average_h_in_m = (Hbar_t +  ...
        SD*pdf('Normal', (x_star(2)/x_star(3)) - Hbar_t, 0, SD)/(1- cdf('Normal', (x_star(2)/x_star(3)) - Hbar_t, 0, SD)));

    Real_world_W_m = x_star(3)*average_h_in_m;

    % Calculate Y_a_t, Y_m_t
    Y_a_t = (Gamma*(Z_lt*x_star(1))^((Mu-1)/Mu) + (1 - Gamma)*(Z_tt*T)^((Mu-1)/Mu))...
        ^(Mu/(Mu - 1));
    Y_m_t = Z_mt * average_h_in_m;

    % Calculate Output per Worker (Total Output in Model since pop = 1)
    Output_perwk = Y_a_t + Y_m_t;

    % Store the results
    endo_t = x_star';
    endo_mat(t,:) = endo_t;

    Real_world_W_m_vec(t) = Real_world_W_m;
    Y_a_t_mat(t) = Y_a_t;
    Y_m_t_mat(t) = Y_m_t;
    Output_perwk_vec(t) = Output_perwk;

end


%% Clean it up (Merging, Renaming)
endo_mat = [endo_mat, Real_world_W_m_vec, Y_a_t_mat, Y_m_t_mat,  Output_perwk_vec];

% Make it into a table
endo_tab = array2table(endo_mat, 'VariableNames', {'L_a', 'W_a', 'W_m', 'W_m_rw', 'Y_a', 'Y_m', 'Out_pwk'});