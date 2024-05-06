cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear

%% Declare parameters values
% Parameters needed are (Z_lt, Z_tt, Z_mt, Gamma, Mu, Hbar_t, SD, P_t, T) %

% Production parameters
Z_lt = 10;
Z_tt = 10;
Z_mt = 80;
Gamma = 0.1;
Mu = 0.5;
T = 1; % Land

%Human Capital parameters
Hbar_t = 0.1;
SD = 1;

% Price
P_t = 2;

%% Try out the function
% Description: % L_a -> x(1), W_a -> x(2), W_m -> x(3)
x0 = [0.5; 20; 2];

x_star = fsolve(@(x)Model2_Function(x, Z_lt, Z_tt, Z_mt, Mu, Gamma, Hbar_t, SD, P_t, T), x0);

%% Try to calculate Land share to see why it's not working as intended.
% Land share = (MPT x T)/Y_a

Y_a = (Gamma*(Z_lt*x_star(1))^((Mu-1)/Mu) + (1 - Gamma)*(Z_tt*T)^((Mu-1)/Mu))...
        ^(Mu/(Mu - 1));
MPT = Z_tt*(1-Gamma)*(Gamma*(Z_lt*x_star(1)/(Z_tt*T))^((Mu -1)/Mu) + (1 - Gamma))^(1/(Mu - 1));

Land_share = (MPT*T)/Y_a;

% It works now. Just as theoritically predicted. 
% Previously, it did not work because land share > Mu. 
%
% With this set of parameters. 
% Mu < Land_share. When increase Z_lt from 10 -> 11, L_a decline. 
