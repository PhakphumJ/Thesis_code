cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
clc

%%% Calibrate to two points (1993 and 2022)
%% Declare parameters values (determined outside of the model)
% Z_m_0, Mu

Z_m_0 = 1; % Normalized to 1 

Mu = 0.5;

g_H_bar =  0.4608; %(growth of mean years of schooling)


%% Exogenous variables 
% Land series
T_series = [6848, 5991];
% Price series
P_series = [1,  1.5309];


%% Set-up the grid of parameters for searching (7 parameters)
% Z_l_0, Z_t_0, H_bar_0
% g_z_m, g_z_l, g_z_t, SD

Z_l_0_list = linspace(0.2,1.2,56);
g_z_l_list = linspace(0.8,1.9,56);

Z_t_0_list = linspace(1.5,6,226);
g_z_t_list = linspace(0.8,2.8,101);

g_z_m_list = linspace(0.8,1.9,56);

H_bar_0_list = linspace(0.1, 4, 40);
SD_list = linspace(0.1, 6, 60);

%% Declare target moments

L_a0 = 0.4239;
L_a1 = 0.3039;

W_a_0_to_W_m_rw_0 = 0.3167;
W_a_1_to_W_m_rw_1 = 0.4178;

Agri_VA_Share_0 = 0.0860;
Agri_VA_Share_1 = 0.0872;

GDP2022toGDPto1993 = 2.6667;

Actual_MM = [L_a0, W_a_0_to_W_m_rw_0, Agri_VA_Share_0, L_a1, W_a_1_to_W_m_rw_1, Agri_VA_Share_1, GDP2022toGDPto1993];

%% Weight Matrix. Give less weight to relative wage since it is less reliable and not the main focus.
Weight = [1.5, 1.25, 1.5, 1.5, 1.25, 1.5, 1.5];

%% Do random serch.
% Draw random combination of parameters. Do it 100,000 times. 
% There are 16,022,476,166,400 possible combination of parameter values. 
n = 5000; %(5,000 for now)

% set seed
rng(2024);

% Matrix to store results.
Search_results_Mat = zeros(n,8);

i = 1;
while i <= n
    % Randomly pick parameter values
    Z_l_0 = randsample(Z_l_0_list,   1);
    g_z_l = randsample(g_z_l_list, 1);
    g_z_m = randsample(g_z_m_list, 1);
    Z_t_0 = randsample(Z_t_0_list, 1);
    g_z_t = randsample(g_z_t_list, 1);
    H_bar_0 = randsample(H_bar_0_list, 1);
    SD = randsample(SD_list, 1);
    
    % Record Parameters
    Param_combi = [Z_l_0; g_z_l; g_z_m; Z_t_0; g_z_t; H_bar_0; SD];
    Search_results_Mat(i,1:7) = Param_combi;

    % Create exogenous variables (Z_mt, Z_lt, Z_tt, H_bar)
    Z_m_series = [Z_m_0, Z_m_0 * (1+g_z_m)];
    Z_l_series = [Z_l_0, Z_l_0 * (1+g_z_l)];
    Z_t_series = [Z_t_0, Z_t_0 * (1+g_z_t)];
    H_bar_series = [H_bar_0, H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (1*7)
    soln_mat = zeros(2,3);

    gdp_mat = zeros(2,1);

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
        Agri_VA_Share_sim = P_t * Y_a_t/(P_t * Y_a_t + Y_m_t);

        % Calculate the norminal gdp
        gdp_sim = P_t * Y_a_t + Y_m_t;

        % Store in the matrix.
        soln_mat(t,:) = [L_a_sim, W_a_to_W_m_rw_sim, Agri_VA_Share_sim];
        gdp_mat(t) = gdp_sim;
    end
        % Convert to 1*6 matrix.
        soln_mat = reshape(soln_mat', [1,6]);

        % Calculate and store the last moment (ratio b/w gdp of the 2 years)
        GDP2022toGDPto1993_sim = gdp_mat(2)/gdp_mat(1);

        soln_mat = [soln_mat, GDP2022toGDPto1993_sim];

        % Calculate Loss (percentage difference squared)
        Loss_mat = ((soln_mat - Actual_MM)./Actual_MM).^2;

        % Summarize by using weighted average
        avg_loss = mean(Weight.*Loss_mat,"all");

        % Store the loss.
        Search_results_Mat(i,8) = avg_loss;

    i = i+1;
end  

%% See what's the best combination of parameters.
% get location
[minloss, index] = min(Search_results_Mat(:,8));

% get the parameter values.
Best_para = Search_results_Mat(index,:);

% Make it into a nice table
Best_para_tab = array2table(Best_para, "VariableNames", {'Z_l_0'; 'g_z_l'; 'g_z_m'; 'Z_t_0'; 'g_z_t'; 'H_bar_0'; 'SD'; 'Loss'});

%% Solve the model with the calibrated parameters.
Best_Z_l_0 = Best_para(1);
Best_g_z_l = Best_para(2);
Best_g_z_m = Best_para(3);
Best_Z_t_0 = Best_para(4);
Best_g_z_t = Best_para(5);
Best_H_bar_0 = Best_para(6);
Best_SD = Best_para(7);

% Create exogenous variables (Z_mt, Z_ll, Z_tt, H_bar)
Z_m_series = [Z_m_0, Z_m_0 * (1+Best_g_z_m)];
Z_l_series = [Best_Z_l_0, Best_Z_l_0 * (1+Best_g_z_l)];
Z_t_series = [Best_Z_t_0, Best_Z_t_0 * (1+Best_g_z_t)];
H_bar_series = [Best_H_bar_0, Best_H_bar_0*(1+g_H_bar)];

% Matrix to store results in each period (2*3)
soln_mat = zeros(2,3);

gdp_mat = zeros(2,1);

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
    x_star = fsolve(@(x)Final_Model_Function(x, Z_lt, Z_tt, Z_mt, Mu, Hbar_t, Best_SD, P_t, T_t), x0);

    % Calculate W_m * avg h in sector m
    h_in_m = (Hbar_t +  ...
    Best_SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/Best_SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/Best_SD, 0, 1)));

    Real_world_W_m = x_star(3)*h_in_m;

    % Calculate Y_a_t, Y_m_t
    Y_a_t = ((Z_lt*x_star(1))^((Mu-1)/Mu) + (Z_tt*T_t)^((Mu-1)/Mu))...
    ^(Mu/(Mu - 1));
    Y_m_t = Z_mt * h_in_m;

    % Calculate the simulated moments
    L_a_sim = x_star(1);
    W_a_to_W_m_rw_sim = x_star(2)/Real_world_W_m;
    Agri_VA_Share_sim = P_t * Y_a_t/(P_t * Y_a_t + Y_m_t);

    % Calculate the norminal gdp
    gdp_sim = P_t * Y_a_t + Y_m_t;

    % Store in the matrix.
    soln_mat(t,:) = [L_a_sim, W_a_to_W_m_rw_sim, Agri_VA_Share_sim];
    gdp_mat(t) = gdp_sim;
end

% Convert to 1*6 matrix.
soln_mat = reshape(soln_mat', [1,6]);

% Calculate and store the last moment (ratio b/w gdp of the 2 years)
GDP2022toGDPto1993_sim = gdp_mat(2)/gdp_mat(1);

soln_mat_calibrated = [soln_mat, GDP2022toGDPto1993_sim];
