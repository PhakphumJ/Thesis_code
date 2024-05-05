cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear

%% Declare parameters values
% Parameters needed are (Z_lt, Z_tt, Z_mt, Gamma, Mu, Hbar_t, SD, P_t, T) %

% Production parameters
Z_lt = 20;
Z_tt = 10;
Z_mt = 50;
Gamma = 0.5;
Mu = 0.7;
T = 10; % Land

%Human Capital parameters
Hbar_t = 2;
SD = 1;

% Price
P_t = 5;

%% Try out the function
% Description: % L_a -> x(1), W_a -> x(2), W_m -> x(3)
x0 = [0.5; 20; 2];

x_star = fsolve(@(x)Model2_Function(x, Z_lt, Z_tt, Z_mt, Mu, Gamma, Hbar_t, SD, P_t, T), x0);