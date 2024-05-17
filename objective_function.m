%% Distance function to be minimized.
% This function would call the fsolve function to solve the model.

function loss = objective_function(params, P_series, T_series, target_vec)
    % Unpack the parameters (from the initial guess of parameters)
    Z_l_0 = params(1);
    g_z_l = params(2);
    g_z_m = params(3);
    Z_t_0 = params(4);
    g_z_t = params(5);
    H_bar_0 = params(6);
    SD = params(7);

    % Create exogenous variables (Z_mt, Z_lt, Z_tt, H_bar)
    Z_m_series = [Z_m_0, Z_m_0 * (1+g_z_m)];
    Z_l_series = [Z_l_0, Z_l_0 * (1+g_z_l)];
    Z_t_series = [Z_t_0, Z_t_0 * (1+g_z_t)];
    H_bar_series = [H_bar_0, H_bar_0*(1+g_H_bar)];

    % Matrix to store results in each period (1*7)
    soln_mat = zeros(2,3);

    gdp_mat = zeros(2,1);

    % Loop over the periods
    for t = 1:2
        % Get the value of exo var.
        P_t = P_series(t);
        T_t = T_series(t);

        Z_lt = Z_l_series(t);
        Z_tt = Z_t_series(t);
        Z_mt = Z_m_series(t);

        Hbar_t = H_bar_series(t);

        % Choose a reasonable starting point of endo var.
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
        Agri_VA_Share_sim = (P_t * Y_a_t)/(P_t * Y_a_t + Y_m_t);

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
    Loss_mat = ((soln_mat - target_vec)./target_vec).^2;

    % Summarize by using weighted average
    avg_loss = mean(Weight.*Loss_mat,"all");
    
    % Return the error
    loss = avg_loss;
end