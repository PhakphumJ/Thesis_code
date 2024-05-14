cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
clc

%%% Calibrate to two points (1993 and 2022)
%% Declare parameters values (determined outside of the model)
% Z_l_0, Z_t_0, g_z_l, g_z_t, g_H_bar, g_P, Mu, P_0

Z_l_0 = 10;
g_z_l = 0.5; %(growth over the 30 years period)

Z_t_0 = 5;
g_z_t = 0.1;

g_H_bar = 0.25; %(growth of mean years of schooling)


Mu = 0.5;



%% Exogenous variables 
% Land series
T_series = [10, 11];

% Price series
P_series = [1 ,0.98];

% Z_l series
Z_l_series = [Z_l_0, Z_l_0*(1+ g_z_l)];

% Z_t series
Z_t_series = [Z_t_0, Z_t_0*(1 + g_z_t)];

%% Set-up the grid of parameters for searching
% Z_m_0, H_bar_0, SD
% g_z_m

Z_m_0_list = linspace(30,130, 101);
g_z_m_list = linspace(0.15,3.5,336);

H_bar_0_list = linspace(50,1000,476);
SD_list = linspace(2,400,200);

%% Declare target moments

L_a0 = 0.5;
L_a1 = 0.3;

W_a_0_to_W_m_rw_0 = 0.3;
W_a_1_to_W_m_rw_1 = 0.5;

Y_a_0_to_Y_0 = 0.15;
Y_a_1_to_Y_1 = 0.10;

Actual_MM = [L_a0, L_a1, W_a_0_to_W_m_rw_0; W_a_1_to_W_m_rw_1, Y_a_0_to_Y_0, Y_a_1_to_Y_1];

%% Do random serch.
% Draw random combination of parameters. Do it 3,000 times.
n = 3000;

% set seed
rng(2024);

% Matrix to store results.
Search_results_Mat = zeros(n,5);

i = 1;
while i <= n
    % Randomly pick parameter values
    Z_m_0 = randsample(Z_m_0_list, 1);
    g_z_m = randsample(g_z_m_list, 1);
    H_bar_0 = randsample(H_bar_0_list, 1);
    SD = randsample(SD_list, 1);
    
    % Record Parameters
    Param_combi = [Z_m_0; g_z_m; H_bar_0; SD];
    Search_results_Mat(i,1:4) = Param_combi;

    % Create exogenous variables (Z_mt, H_bar)
    Z_m_series = [Z_m_0, Z_m_0 * (1+g_z_m)];
    H_bar_series = [H_bar_0, H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each perio  d (2*3)
    soln_mat = zeros(2,3);

    % Solve the model for the two period.
    for t = 1:2
        Z_lt = Z_l_series(t);
        Z_tt = Z_t_series(t);
        Z_mt = Z_m_series(t);

        Hbar_t = H_bar_series(t);
        P_t = P_series(t);
        T_t = T_series(t);
        
        % Choose a reasonable starting point
        x0 = [0.5; Z_lt; Z_mt];

        % Find the equilibrium
        x_star = fsolve(@(x)Final_Model_Function(x, Z_lt, Z_tt, Z_mt, Mu, Hbar_t, SD, P_t, T_t), x0);

        % Calculate W_m * avg h in sector m
        h_in_m = (Hbar_t +  ...
        SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)));

        Real_world_W_m = x_star(3)*h_in_m;

        % Calculate Y_a_t, Y_m_t
        Y_a_t = ((Z_lt*x_star(1))^((Mu-1)/Mu) + (Z_tt*T_t)^((Mu-1)/Mu))...
        ^(Mu/(Mu - 1));
        Y_m_t = Z_mt * h_in_m;

        % Calculate the simulated moments
        L_a_sim = x_star(1);
        W_a_to_W_m_rw_sim = x_star(2)/Real_world_W_m;
        Y_a_to_Y_sim = Y_a_t/(Y_a_t + Y_m_t);

        % Store in the matrix.
        soln_mat(t,:) = [L_a_sim, W_a_to_W_m_rw_sim, Y_a_to_Y_sim];
    end
        %

    i = i+1;
end


