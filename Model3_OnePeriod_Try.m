cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
%% Declare parameters values
% Parameters needed are (Z_lt, Z_tt, Z_mt, Mu, Hbar_t, SD, P_t, T_t) %

% Production parameters
Z_lt = 10;
Z_tt = 5;
Z_mt = 30;
Mu = 0.5;
T_t = 1; % Land

%Human Capital parameters
Hbar_t = 0.1;
SD = 1;

% Price
P_t = 1;

%% Try out the function
% Description: % L_a -> x(1), W_a -> x(2), W_m -> x(3)
x0 = [0.5; 20; 2];

x_star = fsolve(@(x)Final_Model_Function(x, Z_lt, Z_tt, Z_mt, Mu, Hbar_t, SD, P_t, T_t), x0);

%% Try to calculate Land share to see why it's not working as intended.
% Land share = (MPT x T)/Y_a

Y_a = ((Z_lt*x_star(1))^((Mu-1)/Mu) + (Z_tt*T_t)^((Mu-1)/Mu))...
        ^(Mu/(Mu - 1));
MPT = Z_tt^((Mu-1)/Mu)*T_t^(-1/Mu)*((Z_lt*x_star(1))^...
    ((Mu-1)/Mu) + (Z_tt*T_t)^((Mu-1)/Mu))^(1/(Mu - 1));

Land_share = (MPT*T_t)/Y_a;

% It works now. Just as theoritically predicted. 
% Previously, it did not work because land share > Mu. 
%
% With this set of parameters. 
% Mu < Land_share. When increase Z_lt from 10 -> 11, L_a decline. 