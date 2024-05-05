cd 'D:\OneDrive - Central European University\CEU\Thesis\Thesis_code'

clear
%% Declare parameters values (determined outside of the model)
% Z_a_0, g_z_a, alpha, H_bar_0, g_P

Z_a_0 = 5;
alpha = 0.3;
g_z_a = 0.1;

H_bar_0 = 1;
g_P = -0.01;

T = 5; %Land

N = 30; %Num of periods

%% Set-up the grid of parameters for searching
% Z_m_0, g_z_m, g_H_bar, SD, P_0 (6 parameters)

Z_m_0_list = linspace(30,100,2);
g_z_m_list = linspace(0,0.8,2);
g_H_bar_list = linspace(0.01,0.1,2);
SD_list = linspace(0.1,2.5,2);
P_0_list = linspace(0, 1, 2);
H_bar_0_list = linspace(0.1, 10, 2);


%% Declare target moments (only need 6)
% All moments: L_a1, L_a30, W_a_1/W_m_rw_1, W_a_30/W_m_rw_30, Y_a_1/Y_1, Y_a_30/Y_30,
% Out_pwk_1, Out_pwk_30. Choose 5.

L_a1 = 0.5;
L_a30 = 0.3;

W_a_1_to_W_m_rw_1 = 0.3;
W_a_30_to_W_m_rw_30 = 0.5;

Y_a_1_to_Y_1 = 0.15;
Y_a_30_to_Y_30 = 0.10;

Actual_MM = [L_a1, L_a30, W_a_1_to_W_m_rw_1, W_a_30_to_W_m_rw_30, Y_a_1_to_Y_1, Y_a_30_to_Y_30];

%% Calibrate
options = optimoptions(@fsolve, 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, 'FunctionTolerance', 1.0e-05, 'Algorithm', 'levenberg-marquardt');

% Initialize Loss
Loss = 10000;

for Z_m_0 = Z_m_0_list
    for g_z_m = g_z_m_list
        for g_H_bar = g_H_bar_list
            for SD = SD_list
                for P_0 = P_0_list
                    for H_bar_0 = H_bar_0_list
                        % Create Exogenous Variable
                        Z_a = [];
                        Z_a = [Z_a; Z_a_0];
                        Z_m = [];
                        Z_m = [Z_m; Z_m_0];
                        P = [];
                        P = [P; P_0];
                        H_bar = [];
                        H_bar = [H_bar; H_bar_0];

                        for t = 1:N-1
                            z_a_t = Z_a(t);
                            Z_a_t_plus1 = z_a_t*(1 + g_z_a);
                            Z_a = [Z_a; Z_a_t_plus1];
    
                            z_m_t = Z_m(t);
                            Z_m_t_plus1 = z_m_t*(1 + g_z_m);
                            Z_m = [Z_m; Z_m_t_plus1];
    
                            p_t = P(t);
                            p_t_plus1 = p_t*(1 + g_P);
                            P = [P; p_t_plus1];
    
                            H_bar_t = H_bar(t);
                            H_bar_t_plus1 = H_bar_t*(1 + g_H_bar);
                            H_bar = [H_bar; H_bar_t_plus1];    
                        end
    
    
                        % Solve the model
                        try
                            x0 = [0.5; 2; 10]; % Just a starting point for fsolve.

                            endo_mat = zeros(N,3); % To store results from all t.
                            Real_world_W_m_vec = zeros(N,1); %To store results from all t.
                            Y_a_t_mat = zeros(N,1); %To store results from all t.
                            Y_m_t_mat = zeros(N,1); %To store results from all t.
                            Output_perwk_vec = zeros(N,1); %To store results from all t.

                            for t = 1:N
                                Z_at = Z_a(t);
                                Z_mt = Z_m(t);
                                Hbar_t = H_bar(t);
                                P_t = P(t);
                            
                                % Find the equilibrium
                                x_star = fsolve(@(x)Model1_Function(x, Z_at, Z_mt, alpha, Hbar_t, SD, P_t, T), x0, options);
                                
                                % Calculate W_m * avg h in sector m
                                average_h_in_m = (Hbar_t +  ...
                                    SD*pdf('Normal', (x_star(2)/x_star(3)) - Hbar_t, 0, SD)/(1- cdf('Normal', (x_star(2)/x_star(3)) - Hbar_t, 0, SD)));
                            
                                Real_world_W_m = x_star(3)*average_h_in_m;
                            
                                % Calculate Y_a_t, Y_m_t
                                Y_a_t = Z_at*(T^alpha)*x_star(1)^(1-alpha);
                                Y_m_t = Z_mt * average_h_in_m;
                            
                                % Calculate Output per Worker (Total Output in Model since pop = 1)
                                Output_perwk = Y_a_t + Y_m_t;
                            
                                % Store the results
                                endo_t = x_star';
                                endo_mat(t,:) = endo_t;
                            
                                Real_world_W_m_vec(t) = Real_world_W_m;
                                Y_a_t_mat(t) = Y_a_t;
                                Y_m_t_mat(t) = Y_m_t;
                                Output_perwk_vec(t) = Output_perwk;
                            end
                            
                            % Compile the relavant endogenous
                            % Description: L_a -> x(1), W_a -> x(2), W_m -> x(3)
                            Generated_MM = [endo_mat(1,1), endo_mat(30,1), ...
                            Real_world_W_m_vec(1), Real_world_W_m_vec(30), Y_a_t_mat(1)/Output_perwk(1),...
                            Y_a_t_mat(30)/Output_perwk(30)];

                            % Calculate the loss (percentage difference squared)
                            Loss_vector = (log(Generated_MM) - log(Actual_MM)).^2;

                            new_loss = mean(Loss_vector); % Hence, equal weight for each moment.
                        catch
                            fprintf('Found complex root, skip \n');
                        end
                        
                        if new_loss < Loss
                            Cor_Para = [Z_m_0, g_z_m, g_H_bar, SD, P_0, H_bar_0];
                            Loss = new_loss;
                        end
                    end
                end
            end
        end
    end
end

