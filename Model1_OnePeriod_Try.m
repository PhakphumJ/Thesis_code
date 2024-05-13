cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
%% Declare parameters values
% Parameters needed are (Z_at, Z_mt, alpha, Hbar_t, SD, P_t, T) %

% Production parameters
Z_at = 10;
Z_mt = 50;
alpha = 0.3;
T = 1000; % Land

%Human Capital parameters
Hbar_t = 5;
SD = 1;

% Price
P_t = 5;

%% Try out the function
% Description: % L_a -> x(1), W_a -> x(2), W_m -> x(3)
x0 = [0.5; 20; 2];

x_star = fsolve(@(x)Model1_Function(x, Z_at, Z_mt, alpha, Hbar_t, SD, P_t, T), x0);

