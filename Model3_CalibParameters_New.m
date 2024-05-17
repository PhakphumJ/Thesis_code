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
T_series = [6.8489, 5.9502];

% Price series
P_series = [0.9247,  1.4156];


%% Set-up the grid of parameters for searching (7 parameters)
% Z_l_0, Z_t_0, H_bar_0
% g_z_m, g_z_l, g_z_t, SD

Z_l_0_list = linspace(0.1,1.3,61);
g_z_l_list = linspace(0.4,2.1,106);

Z_t_0_list = linspace(0.1,1.3,61);
g_z_t_list = linspace(0.4,1.9,86);

g_z_m_list = linspace(0.4,1.9,86);

H_bar_0_list = linspace(0.1, 5, 50);
SD_list = logspace(log10(0.05), log10(15), 175);

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
Weight = [1, 1, 1, 1, 1, 1, 1];

%% Do random serch.
% Draw random combination of parameters. Do it 500,000 times. 
% There are 6,543,049,590,000 possible combination of parameter values. 
n = 10000; %(5,000 for now)

% set seed
rng(2024);

% Matrix to store results.
Search_results_Mat = zeros(n,8);

i = 1;
while i <= n
    % Randomly pick parameter values
    Z_l_0 = randsample(Z_l_0_list, 1);
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
        soln_mat = reshape(soln_mat, [1,6]);

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
soln_mat = reshape(soln_mat, [1,6]);

% Calculate and store the last moment (ratio b/w gdp of the 2 years)
GDP2022toGDPto1993_sim = gdp_mat(2)/gdp_mat(1);

soln_mat_calibrated = [soln_mat, GDP2022toGDPto1993_sim];


%% Counterfactual analysis (change g_z_l, and g_H_bar)
% Case 1: Only increase g_H_bar
% Case 2: Only increse g_z_l
% Case 3: Increase both

new_g_z_l = g_z_l*1.20;
new_g_H_bar = Best_g_H_bar*1.20;

% Only need to solve for the last period.
% Create exogenous variables (Z_l_1, Z_t_1, Z_m_1, H_bar_1, P_1, T_1)

new_Z_l_1 = Z_l_0*(1 + new_g_z_l);
old_Z_l_1 = Z_l_0*(1 + g_z_l);
Z_t_1 = Z_t_series(2);
Z_m_1 = Z_m_series(2);
new_H_bar_1 = Best_H_bar_0*(1 + new_g_H_bar);
old_H_bar_1 = Best_H_bar_0*(1 + Best_g_H_bar);
P_1 = P_series(2);
T_1 = T_series(2);

% Solve the model for the counterfactual parameters.
% Find the equilibrium
%% Case 1 (change only g_H_bar)
x_star_counter_case1 = fsolve(@(x)Final_Model_Function(x, old_Z_l_1, Z_t_1, Z_m_1, Mu, new_H_bar_1, SD, P_1, T_1), x0);

% Calculate W_m * avg h in sector m
h_in_m = (new_H_bar_1 +  ...
SD*pdf('Normal', ((x_star_counter_case1(2)/x_star_counter_case1(3)) - new_H_bar_1)/SD, 0, 1)/(1- cdf('Normal', ((x_star_counter_case1(2)/x_star_counter_case1(3)) - new_H_bar_1)/SD, 0, 1)));

Real_world_W_m = x_star_counter_case1(3)*h_in_m;

% Calculate Y_a_t, Y_m_t
Y_a_t = ((old_Z_l_1*x_star_counter_case1(1))^((Mu-1)/Mu) + (Z_t_1*T_1)^((Mu-1)/Mu))...
^(Mu/(Mu - 1));
Y_m_t = Z_m_1 * h_in_m;

% Calculate the simulated moments
L_a_sim = x_star_counter_case1(1);
W_a_to_W_m_rw_sim = x_star_counter_case1(2)/Real_world_W_m;
Agri_VA_Share_sim = Y_a_t/(Y_a_t + Y_m_t);

% Store in the matrix.
soln_mat_counter_1 = [L_a_sim, W_a_to_W_m_rw_sim, Agri_VA_Share_sim];

%% Case 2 (change only g_z_l)
x_star_counter_case2 = fsolve(@(x)Final_Model_Function(x, new_Z_l_1, Z_t_1, Z_m_1, Mu, old_H_bar_1, SD, P_1, T_1), x0);

% Calculate W_m * avg h in sector m
h_in_m = (old_H_bar_1 +  ...
SD*pdf('Normal', ((x_star_counter_case2(2)/x_star_counter_case2(3)) - old_H_bar_1)/SD, 0, 1)/(1- cdf('Normal', ((x_star_counter_case2(2)/x_star_counter_case2(3)) - old_H_bar_1)/SD, 0, 1)));

Real_world_W_m = x_star_counter_case2(3)*h_in_m;

% Calculate Y_a_t, Y_m_t
Y_a_t = ((new_Z_l_1*x_star_counter_case2(1))^((Mu-1)/Mu) + (Z_t_1*T_1)^((Mu-1)/Mu))...
^(Mu/(Mu - 1));
Y_m_t = Z_m_1 * h_in_m;

% Calculate the simulated moments
L_a_sim = x_star_counter_case2(1);
W_a_to_W_m_rw_sim = x_star_counter_case2(2)/Real_world_W_m;
Agri_VA_Share_sim = Y_a_t/(Y_a_t + Y_m_t);

% Store in the matrix.
soln_mat_counter_2 = [L_a_sim, W_a_to_W_m_rw_sim, Agri_VA_Share_sim];

%% Case 3 (change both of the two growth rates)
x_star_counter_case3 = fsolve(@(x)Final_Model_Function(x, new_Z_l_1, Z_t_1, Z_m_1, Mu, new_H_bar_1, SD, P_1, T_1), x0);

% Calculate W_m * avg h in sector m
h_in_m = (new_H_bar_1 +  ...
SD*pdf('Normal', ((x_star_counter_case3(2)/x_star_counter_case3(3)) - new_H_bar_1)/SD, 0, 1)/(1- cdf('Normal', ((x_star_counter_case3(2)/x_star_counter_case3(3)) - new_H_bar_1)/SD, 0, 1)));

Real_world_W_m = x_star_counter_case3(3)*h_in_m;

% Calculate Y_a_t, Y_m_t
Y_a_t = ((new_Z_l_1*x_star_counter_case3(1))^((Mu-1)/Mu) + (Z_t_1*T_1)^((Mu-1)/Mu))...
^(Mu/(Mu - 1));
Y_m_t = Z_m_1 * h_in_m;

% Calculate the simulated moments
L_a_sim = x_star_counter_case3(1);
W_a_to_W_m_rw_sim = x_star_counter_case3(2)/Real_world_W_m;
Agri_VA_Share_sim = Y_a_t/(Y_a_t + Y_m_t);

% Store in the matrix.
soln_mat_counter_3 = [L_a_sim, W_a_to_W_m_rw_sim, Agri_VA_Share_sim];

%% Compile the results (Actual data, Baseline, Counter Factual Case 1, 2, and 3)

Compiled_results_mat = [Actual_MM(2,:); soln_mat_calibrated(2,:); soln_mat_counter_1; soln_mat_counter_2; soln_mat_counter_3];

Compiled_results_tab = array2table(Compiled_results_mat, "VariableNames", {'Agricultural Employment Share', 'Relative Wages', 'Agricultural Value-added Share'});

Compiled_results_tab.Properties.RowNames = {'Data', 'Baseline', 'Increased growth rate of H_bar', 'Increased growth rate of Z_l', 'Incresed growth rates of H_bar and Z_l'};
