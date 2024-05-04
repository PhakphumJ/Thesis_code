cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
%% Declare parameters values
% Production parameters
Z_a_0 = 5;
Z_m_0 = 40;
alpha = 0.3;
g_z_a = 0.07;
g_z_m = 0.05;

% Human capital parameters
H_bar_0 = 1;
g_H_bar = 0.05;
SD = 1.5;

% Utility parameters
Beta = 0.98;

% Price
P_0 = 5;
g_P = 0.01;


N = 30;
%% Create Exogeneous Variables (Z, p, H_bar)
Z_a = [];
Z_a = [Z_a; Z_a_0];
for t = 1:N-1
   z_a_t = Z_a(t);
   Z_a_t_plus1 = z_a_t*(1 + g_z_a);
   Z_a = [Z_a; Z_a_t_plus1];
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

T = 5; %Land
%% Solve the model for each t
% Description: % L_a -> x(1), Y_a -> x(2), Y_m -> x(3), W_a -> x(4), 
% W_m -> x(5)

x0 = [0.5; 20; 2; 10; 30]; % Just a starting point for fsolve.

endo_mat = zeros(N,5); % To store results from all t.
Real_world_W_m_vec = zeros(N,1); %To store results from all t.

for t = 1:N
    Z_at = Z_a(t);
    Z_mt = Z_m(t);
    Hbar_t = H_bar(t);
    P_t = P(t);

    % Find the equilibrium
    x_star = fsolve(@(x)Model1_Function(x, Z_at, Z_mt, alpha, Hbar_t, SD, P_t, T), x0);
    
    % Calculate W_m * avg h in sector m
    Real_world_W_m = x_star(5)*(Hbar_t +  ...
        SD*pdf('Normal', (x_star(4)/x_star(5)) - Hbar_t, 0, SD)/(1- cdf('Normal', (x_star(4)/x_star(5)) - Hbar_t, 0, SD)));

    % Store the results
    endo_t = x_star';
    endo_mat(t,:) = endo_t;

    Real_world_W_m_vec(t) = Real_world_W_m;

end


%% Clean it up (Merging, Renaming)
endo_mat = [endo_mat, Real_world_W_m_vec];

% Make it into a table
endo_tab = array2table(endo_mat, 'VariableNames', {'L_a', 'Y_a', 'Y_m', 'W_a', 'W_m', 'W_m_rw'});