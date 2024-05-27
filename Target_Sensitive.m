% This is for examining the sensitity of target moments with respect to
% parameters. 

%% First let's calibrate the model.
cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
clc

%%% Calibrate to two points (1993 and 2022)

%% Declare parameters values (determined outside of the model)
% Z_m_0, Mu
global Z_m_0
Z_m_0 = 1; % Normalized to 1 

global Mu
Mu = 0.5;

global g_H_bar
g_H_bar =  0.4608; %(growth of mean years of schooling)

%% Exogenous variables 
% Land series
T_series = [6848, 5991];

% Price series
P_series = [1,  1.5309];

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
global Weight
Weight = [1.6, 1, 1.6, 1.6, 1, 1.6, 1.6];


%% Optimize the parameters
% Initial guess for the parameters

initial_params = [0.1, 1, 1, 1, 1, 0.5, 0.5];

% Define lower and upper bounds for the parameters
lb = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05];
ub = [2, 7, 7, 2, 7, 10, 15];  

% use simplex method to select the parameters.
[optimized_params, fval] = fminsearchbnd(@(params) objective_function(params, P_series, T_series, Actual_MM), initial_params, lb, ub);

% Display the optimized parameters
disp('Optimized Parameters:');
disp(optimized_params);

% Make it into a nice table
Best_para_tab = array2table(optimized_params', "RowNames", {'Z_l_2001'; 'g_Z_l'; 'g_Z_m'; 'Z_T_2001'; 'g_Z_T'; 'H_bar_2001'; 'Sigma'});
new_order = [1,2,4,5,3,6,7];
Best_para_tab = Best_para_tab(new_order,:);
Best_para_tab.Properties.VariableNames = {'Values'};

% Export to CSV
writetable(Best_para_tab,'Calibrated_Parameters.csv', 'WriteRowNames',true);

%% Try out the obtained parameters.

% Extract the parameter values.
Best_Z_l_0 = optimized_params(1);
Best_g_z_l = optimized_params(2);
Best_g_z_m = optimized_params(3);
Best_Z_t_0 = optimized_params(4);
Best_g_z_t = optimized_params(5);
Best_H_bar_0 = optimized_params(6);
Best_SD = optimized_params(7);

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

% Calculate land share.
condition = (Z_tt*T_t)^((Mu-1)/Mu)/((Z_lt*x_star(1))^((Mu-1)/Mu) + (Z_tt*T_t)^((Mu-1)/Mu));

% Compile into a table (Actual and Simulated)
Target_compare_tab = array2table([Actual_MM ; soln_mat_calibrated]', "VariableNames", {'Actual', 'Simulated'});

% Add row names
Target_compare_tab.Properties.RowNames = {'Agricultural Employment Share 2001'; 'Relative Wage 2001'; 'Agricultural Value-added Share 2001'; 'Agricultural Employment Share 2022'; 'Relative Wage 2022'; 'Agricultural Value-added Share 2022'; 'GDP 2022/ GDP 2001'};

% Change row order
newOrder = [1, 4, 2, 5, 3, 6, 7]; 
Target_compare_tab = Target_compare_tab(newOrder,:);

%% Sensitivity 
% for parameter with pp as unit. Examine the change in parameter values by
% +- 0.01 pp, +-0.02 pp, +-0.03 pp, +- and 0.04 pp
% for the rest, change the value of parameters by:
% +- 0.05, +-0.10, +-0.15, +- and 0.20

%% Let's start with changing sigma.

% Create the list of parameter being examined.
sigma_list = [Best_SD - 0.20, Best_SD - 0.15, Best_SD - 0.10, Best_SD - 0.05, ...
    Best_SD, Best_SD + 0.05, Best_SD + 0.10, Best_SD + 0.15, Best_SD + 0.20];

% Create a matrix to store results. 7 x 9
sen_sigma_result = zeros(7,9);


% Solve the model and compute the target moments.


for sigma = 1:length(sigma_list)
    SD = sigma_list(sigma);

    % Create exogenous variables. 
    Z_m_series = [Z_m_0, Z_m_0 * (1+Best_g_z_m)];
    Z_l_series = [Best_Z_l_0, Best_Z_l_0 * (1+Best_g_z_l)];
    Z_t_series = [Best_Z_t_0, Best_Z_t_0 * (1+Best_g_z_t)];
    H_bar_series = [Best_H_bar_0, Best_H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (2*3)
    soln_mat = zeros(2,3);
    
    gdp_mat = zeros(2,1);

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

    % Compile the coumputed target moments for this value of parameter.
    sen_sigma_result(:, sigma) = soln_mat;
end

% Compute the % change from the baseline case.
Baseline_case = sen_sigma_result(:,5);

Change_sen_sigma_result = (sen_sigma_result - Baseline_case) *100 ./Baseline_case;


% Make it into a nice table. 
Change_sen_sigma_result_tab = array2table(Change_sen_sigma_result, "RowNames", {'Agricultural Employment Share 2001'; 'Relative Wage 2001'; 'Agricultural Value-added Share 2001'; 'Agricultural Employment Share 2022'; 'Relative Wage 2022'; 'Agricultural Value-added Share 2022'; 'GDP 2022/ GDP 2001'});
Change_sen_sigma_result_tab.Properties.VariableNames = {'-0.20', '-0.15', '-0.10', '0.05', 'Cailibrated Value', '+0.05', '+0.10', '+0.15', '+0.20'};

%% Next, h_bar_0
% Create the list of parameter being examined.
H_bar_0_list = [Best_H_bar_0 - 0.20, Best_H_bar_0 - 0.15, Best_H_bar_0 - 0.10, Best_H_bar_0 - 0.05, ...
    Best_H_bar_0, Best_H_bar_0 + 0.05, Best_H_bar_0 + 0.10, Best_H_bar_0 + 0.15, Best_H_bar_0 + 0.20];

% Create a matrix to store results. 7 x 9
sen_H_bar_0_result = zeros(7,9);


% Solve the model and compute the target moments.


for H_bar_0 = 1:length(H_bar_0_list)
    h_bar_0 = H_bar_0_list(H_bar_0);

    % Create exogenous variables. 
    Z_m_series = [Z_m_0, Z_m_0 * (1+Best_g_z_m)];
    Z_l_series = [Best_Z_l_0, Best_Z_l_0 * (1+Best_g_z_l)];
    Z_t_series = [Best_Z_t_0, Best_Z_t_0 * (1+Best_g_z_t)];
    H_bar_series = [h_bar_0, h_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (2*3)
    soln_mat = zeros(2,3);
    
    gdp_mat = zeros(2,1);

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
        Best_SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/Best_SD, 0, 1)));
    
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

    % Compile the coumputed target moments for this value of parameter.
    sen_H_bar_0_result(:, H_bar_0) = soln_mat;
end

% Compute the % change from the baseline case.
Baseline_case = sen_H_bar_0_result(:,5);
Change_sen_H_bar_0_result = (sen_H_bar_0_result - Baseline_case) *100 ./Baseline_case;

% Make it into a nice table. 
Change_sen_H_bar_0_result_tab = array2table(Change_sen_H_bar_0_result, "RowNames", {'Agricultural Employment Share 2001'; 'Relative Wage 2001'; 'Agricultural Value-added Share 2001'; 'Agricultural Employment Share 2022'; 'Relative Wage 2022'; 'Agricultural Value-added Share 2022'; 'GDP 2022/ GDP 2001'});
Change_sen_H_bar_0_result_tab.Properties.VariableNames = {'-0.20', '-0.15', '-0.10', '0.05', 'Cailibrated Value', '+0.05', '+0.10', '+0.15', '+0.20'};

%% Next, Z_L_0
% Create the list of parameter being examined.
Z_L_0_list = [Best_Z_l_0 - 0.20, Best_Z_l_0 - 0.15, Best_Z_l_0 - 0.10, Best_Z_l_0 - 0.05, ...
    Best_Z_l_0, Best_Z_l_0 + 0.05, Best_Z_l_0 + 0.10, Best_Z_l_0 + 0.15, Best_Z_l_0 + 0.20];

% Create a matrix to store results. 7 x 9
sen_Z_L_0_result = zeros(7,9);


% Solve the model and compute the target moments.


for Z_L_0 = 1:length(Z_L_0_list)
    Z_l_0 = Z_L_0_list(Z_L_0);

    % Create exogenous variables. 
    Z_m_series = [Z_m_0, Z_m_0 * (1+Best_g_z_m)];
    Z_l_series = [Z_l_0, Z_l_0 * (1+Best_g_z_l)];
    Z_t_series = [Best_Z_t_0, Best_Z_t_0 * (1+Best_g_z_t)];
    H_bar_series = [Best_H_bar_0, Best_H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (2*3)
    soln_mat = zeros(2,3);
    
    gdp_mat = zeros(2,1);

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
        Best_SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/Best_SD, 0, 1)));
    
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

    % Compile the coumputed target moments for this value of parameter.
    sen_Z_L_0_result(:, Z_L_0) = soln_mat;
end

% Compute the % change from the baseline case.
Baseline_case = sen_Z_L_0_result(:,5);
Change_sen_Z_L_0_result = (sen_Z_L_0_result - Baseline_case) *100 ./Baseline_case;

% Make it into a nice table. 
Change_sen_Z_L_0_result_tab = array2table(Change_sen_Z_L_0_result, "RowNames", {'Agricultural Employment Share 2001'; 'Relative Wage 2001'; 'Agricultural Value-added Share 2001'; 'Agricultural Employment Share 2022'; 'Relative Wage 2022'; 'Agricultural Value-added Share 2022'; 'GDP 2022/ GDP 2001'});
Change_sen_Z_L_0_result_tab.Properties.VariableNames = {'-0.20', '-0.15', '-0.10', '0.05', 'Cailibrated Value', '+0.05', '+0.10', '+0.15', '+0.20'};

%% Next, Z_T_0
% Create the list of parameter being examined.
Z_T_0_list = [Best_Z_t_0 - 0.20, Best_Z_t_0 - 0.15, Best_Z_t_0 - 0.10, Best_Z_t_0 - 0.05, ...
    Best_Z_t_0, Best_Z_t_0 + 0.05, Best_Z_t_0 + 0.10, Best_Z_t_0 + 0.15, Best_Z_t_0 + 0.20];

% Create a matrix to store results. 7 x 9
sen_Z_T_0_result = zeros(7,9);


% Solve the model and compute the target moments.


for Z_T_0 = 1:length(Z_T_0_list)
    Z_t_0 = Z_T_0_list(Z_T_0);

    % Create exogenous variables. 
    Z_m_series = [Z_m_0, Z_m_0 * (1+Best_g_z_m)];
    Z_l_series = [Best_Z_l_0, Best_Z_l_0 * (1+Best_g_z_l)];
    Z_t_series = [Z_t_0, Z_t_0 * (1+Best_g_z_t)];
    H_bar_series = [Best_H_bar_0, Best_H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (2*3)
    soln_mat = zeros(2,3);
    
    gdp_mat = zeros(2,1);

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
        Best_SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/Best_SD, 0, 1)));
    
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

    % Compile the coumputed target moments for this value of parameter.
    sen_Z_T_0_result(:, Z_T_0) = soln_mat;
end

% Compute the % change from the baseline case.
Baseline_case = sen_Z_T_0_result(:,5);
Change_sen_Z_T_0_result = (sen_Z_T_0_result - Baseline_case) *100 ./Baseline_case;

% Make it into a nice table. 
Change_sen_Z_T_0_result_tab = array2table(Change_sen_Z_T_0_result, "RowNames", {'Agricultural Employment Share 2001'; 'Relative Wage 2001'; 'Agricultural Value-added Share 2001'; 'Agricultural Employment Share 2022'; 'Relative Wage 2022'; 'Agricultural Value-added Share 2022'; 'GDP 2022/ GDP 2001'});
Change_sen_Z_T_0_result_tab.Properties.VariableNames = {'-0.20', '-0.15', '-0.10', '0.05', 'Cailibrated Value', '+0.05', '+0.10', '+0.15', '+0.20'};

%% Next, G_Z_L
% Create the list of parameter being examined.
G_Z_L_list = [Best_g_z_l - 0.04, Best_g_z_l - 0.03, Best_g_z_l - 0.02, Best_g_z_l - 0.01, ...
    Best_g_z_l, Best_g_z_l + 0.01, Best_g_z_l + 0.02, Best_g_z_l + 0.03, Best_g_z_l + 0.04];

% Create a matrix to store results. 7 x 9
sen_G_Z_L_result = zeros(7,9);


% Solve the model and compute the target moments.


for G_Z_L = 1:length(G_Z_L_list)
    G_Z_l = G_Z_L_list(Z_T_0);

    % Create exogenous variables. 
    Z_m_series = [Z_m_0, Z_m_0 * (1+Best_g_z_m)];
    Z_l_series = [Best_Z_l_0, Best_Z_l_0 * (1+G_Z_l)];
    Z_t_series = [Best_Z_t_0, Best_Z_t_0 * (1+Best_g_z_t)];
    H_bar_series = [Best_H_bar_0, Best_H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (2*3)
    soln_mat = zeros(2,3);
    
    gdp_mat = zeros(2,1);

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
        Best_SD*pdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/SD, 0, 1)/(1- cdf('Normal', ((x_star(2)/x_star(3)) - Hbar_t)/Best_SD, 0, 1)));
    
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

    % Compile the coumputed target moments for this value of parameter.
    sen_G_Z_L_result(:, G_Z_L) = soln_mat;
end

% Compute the % change from the baseline case.
Baseline_case = sen_G_Z_L_result(:,5);
Change_sen_G_Z_L_result = (sen_G_Z_L_result - Baseline_case) *100 ./Baseline_case;

% Make it into a nice table. 
Change_sen_G_Z_L_result_tab = array2table(Change_sen_G_Z_L_result, "RowNames", {'Agricultural Employment Share 2001'; 'Relative Wage 2001'; 'Agricultural Value-added Share 2001'; 'Agricultural Employment Share 2022'; 'Relative Wage 2022'; 'Agricultural Value-added Share 2022'; 'GDP 2022/ GDP 2001'});
Change_sen_G_Z_L_result_tab.Properties.VariableNames = {'-0.04 pp', '-0.03 pp', '-0.02 pp', '-0.01 pp', 'Cailibrated Value', '+0.01 pp', '+0.02 pp', '+0.03 pp', '+0.04 pp'};
