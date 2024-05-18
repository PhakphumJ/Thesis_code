cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
clc

%%% Calibrate to two points (1993 and 2022)

%% Declare parameters values (determined outside of the model)
% Z_m_0, Mu
global Z_m_0
Z_m_0 = 1; % Normalized to 1 

global Mu
Mu = 0.7;

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

initial_params = [0.1, 1, 1, 0.1, 1, 0.5, 0.5];

% Define lower and upper bounds for the parameters
lb = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05];
ub = [1, 5, 5, 5, 1, 10, 15];  

% use simplex method to select the parameters.
[optimized_params, fval] = fminsearchbnd(@(params) objective_function(params, P_series, T_series, Actual_MM), initial_params, lb, ub);

% Display the optimized parameters
disp('Optimized Parameters:');
disp(optimized_params);

% Make it into a nice table
Best_para_tab = array2table(optimized_params, "VariableNames", {'Z_l_0'; 'g_z_l'; 'g_z_m'; 'Z_t_0'; 'g_z_t'; 'H_bar_0'; 'SD'});

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

%% Couter factual analysis 
% (increase z_l by : [5%, 10%, 15%, 20%]
% (increase H_bar by : [5%, 10%, 15%, 20%]
% Hence 5 * 5 cases.

% Specify the counter factual parameters.

Base_Z_l_1 = Z_l_series(2);
Base_H_bar_1 = H_bar_series(2);

Z_l_list = [Base_Z_l_1, Base_Z_l_1*1.05, Base_Z_l_1*1.1, Base_Z_l_1*1.15, Base_Z_l_1*1.20];

H_bar_list = [Base_H_bar_1, Base_H_bar_1*1.05, Base_H_bar_1*1., Base_H_bar_1*1.15, Base_H_bar_1*1.20];


% Matrix to store the results. Only solve for the last period.
Counter_results_mat_LA = zeros(5,5);
Counter_results_mat_VA = zeros(5,5);
Counter_results_mat_GDP = zeros(5,5);

for i = 1:5
    for j = 1:5
       Z_lt = Z_l_list(i);
       Hbar_t = H_bar_list(j);


       % Create the exogenous
       Z_tt = Z_t_series(2);
       Z_mt = Z_m_series(2);
       P_t = P_series(2);
       T_t = T_series(2);

       % Solve the model.
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
        Agri_VA_Share_sim = P_t * Y_a_t/(P_t * Y_a_t + Y_m_t);

        % Calculate the norminal gdp
        gdp_sim = P_t * Y_a_t + Y_m_t;

        % Store in the matrix.
        Counter_results_mat_LA(i,j) = L_a_sim;
        Counter_results_mat_VA(i,j) = Agri_VA_Share_sim;
        Counter_results_mat_GDP(i,j) = gdp_sim;
    end
end

% Create Change in GDP mat.
BaselineGDP_mat = ones(1,1).*Counter_results_mat_GDP(1,1);
Change_GDP_mat = (Counter_results_mat_GDP - BaselineGDP_mat)*100./BaselineGDP_mat;

