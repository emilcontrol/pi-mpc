% This script is to simulate and generate data for the case
% With a slight difference between the model used in the 
% MPC and the simulation model.

clear all
[ad,bd,A,B_u,B_d,R,Ts,N,Ndist,K,Ti,Kp,Ki,K_fb,distVec,distVec_ctrl,x0] = readParamsToWorkspace();

% PI controllers
PI_smaller = cPI();
PI_normal = cPI();
PI_larger = cPI();

% Get Q and Qf
[Q, Qf] = calculateQdiscForPI(ad,bd,Ki,Kp,R);

controlHorizon = 60;
predHorizon = controlHorizon;

% Reference only for plotting. Always zero.
r = [zeros(N, 1), zeros(N, 1)];

% Initial states
% Minus since we here have y and not e
y_pid_smaller = -1*x0(2);
x_mpc_smaller = x0;
y_pid_normal = -1*x0(2);
x_mpc_normal = x0;
y_pid_larger = -1*x0(2);
x_mpc_larger = x0;

% History
u_history_pid_smaller = zeros(N, 1);
y_history_pid_smaller = zeros(N, 1);
u_history_mpc_smaller = zeros(N, 1);
x_history_mpc_smaller = zeros(N, 2);

u_history_pid_normal = zeros(N, 1);
y_history_pid_normal = zeros(N, 1);
u_history_mpc_normal = zeros(N, 1);
x_history_mpc_normal = zeros(N, 2);

u_history_pid_larger = zeros(N, 1);
y_history_pid_larger = zeros(N, 1);
u_history_mpc_larger = zeros(N, 1);
x_history_mpc_larger = zeros(N, 2);

% Optimization matrices H and fnx. G and F we need for constraint matrices
[H,fnx_u,fnx_d, G_u, G_d,F] = CalculateQPMtx(A,B_u,B_d,Q,R,Qf,controlHorizon,predHorizon);

% Hard constraints only for MPC. Note that this is constraints on e=-y. Not y.
Xbounds = [-inf  inf ; % x1min, x1max
           -inf  inf]; % x2min, x2max
Ubounds = [-inf  inf]; % umin, umax

% Constraint matrices.
[Cbar, cnx, Fext, Gdext] = CalculateConstraintMtx(A, B_u, ...
                            Xbounds, Ubounds, controlHorizon, predHorizon, G_u, G_d, F);

% Suppress outputs from quadprog
OptOptions = optimoptions('quadprog', 'Display', 'off','OptimalityTolerance',eps, ...
    'MaxIterations',1000);
warning('off', 'all');

% Shift A for simulation a little bit to make sure that
% model does not match reality. A 5% smaller
A_sim = A;
A_sim(2,2) = A_sim(2,2)*1;


% Shift B for simulation a little bit to make sure that
% model does not match reality. B 10% larger
B_u_sim_smaller = B_u;
B_u_sim_smaller(2) = B_u_sim_smaller(2)*0.50;

B_u_sim_normal = B_u;
B_u_sim_normal(2) = B_u_sim_normal(2)*1;

B_u_sim_larger = B_u;
B_u_sim_larger(2) = B_u_sim_larger(2)*1.50;

% Simulating
for k = 1:N
    %Measurable disturbances
    D = distVec_ctrl(k:k+predHorizon-1);

    % control signal for PID
    [u_pid_smaller, PI_smaller] = PI_smaller.PI(r(k,2),y_pid_smaller,Kp,Ki);
    [u_pid_normal, PI_normal] = PI_normal.PI(r(k,2),y_pid_normal,Kp,Ki);
    [u_pid_larger, PI_larger] = PI_larger.PI(r(k,2),y_pid_larger,Kp,Ki);

    % Solve QP to find optimal U
    % min U'*H*U/2 + x'*fnx*U subject to Cbar*U <= cnx - Fext*x
    tic
    if k == 1
        U_mpc_old_smaller = zeros(controlHorizon,1);
        U_mpc_old_normal = zeros(controlHorizon,1);
        U_mpc_old_larger = zeros(controlHorizon,1);
    else
        U_mpc_old_smaller = circshift(U_mpc_smaller,-1);
        U_mpc_old_normal = circshift(U_mpc_normal,-1);
        U_mpc_old_larger = circshift(U_mpc_larger,-1);
    end
    U_mpc_smaller = SolveMPC(x_mpc_smaller,H,fnx_u,fnx_d, Cbar,cnx,Fext,Gdext,D,U_mpc_old_smaller, OptOptions);
    U_mpc_normal = SolveMPC(x_mpc_normal,H,fnx_u,fnx_d, Cbar,cnx,Fext,Gdext,D,U_mpc_old_normal, OptOptions);
    U_mpc_larger = SolveMPC(x_mpc_larger,H,fnx_u,fnx_d, Cbar,cnx,Fext,Gdext,D,U_mpc_old_larger, OptOptions);
    
    OptTime = toc;

    disp(['Iteration: ', num2str(k), '. Optimization time: ', ...
        num2str(round(OptTime*1000, 1 - floor(log10(abs(OptTime*1000))))), ' ms']);

    % Save simulation states and control inputs
    u_history_mpc_smaller(k) = U_mpc_smaller(1,:);
    x_history_mpc_smaller(k, :) = x_mpc_smaller;
    u_history_pid_smaller(k) = u_pid_smaller;
    y_history_pid_smaller(k) = y_pid_smaller;

    u_history_mpc_normal(k) = U_mpc_normal(1,:);
    x_history_mpc_normal(k, :) = x_mpc_normal;
    u_history_pid_normal(k) = u_pid_normal;
    y_history_pid_normal(k) = y_pid_normal;

    u_history_mpc_larger(k) = U_mpc_larger(1,:);
    x_history_mpc_larger(k, :) = x_mpc_larger;
    u_history_pid_larger(k) = u_pid_larger;
    y_history_pid_larger(k) = y_pid_larger;

    % Disturbance vector
    d=distVec(k);
    
    % Simulate the plant for one step for PID. -1 before B_d since it is created for error=-y.
    % Also -1 before B_u since B_u is created for x (error), which differs
    % on -1
    y_pid_smaller = A_sim(2,2) * y_pid_smaller - B_u_sim_smaller(2) * u_pid_smaller - B_d(2)*d; 
    y_pid_normal = A_sim(2,2) * y_pid_normal - B_u_sim_normal(2) * u_pid_normal - B_d(2)*d; 
    y_pid_larger = A_sim(2,2) * y_pid_larger - B_u_sim_larger(2) * u_pid_larger - B_d(2)*d; 

    % Simulate the plant for one step for MPC
    x_mpc_smaller = A_sim * x_mpc_smaller + B_u_sim_smaller * U_mpc_smaller(1,:) + B_d*d; 
    x_mpc_normal = A_sim * x_mpc_normal + B_u_sim_normal * U_mpc_normal(1,:) + B_d*d; 
    x_mpc_larger = A_sim * x_mpc_larger + B_u_sim_larger * U_mpc_larger(1,:) + B_d*d; 
end

% Translate MPC x=e back to y=r-e
y_history_mpc_smaller = -1*x_history_mpc_smaller(:, 2);
y_history_mpc_normal = -1*x_history_mpc_normal(:, 2);
y_history_mpc_larger = -1*x_history_mpc_larger(:, 2);

% PID
%e_history_pid_smaller = r(:,1) - y_history_pid_smaller;
%e_history_pid_normal = r(:,1) - y_history_pid_normal;
%e_history_pid_larger = r(:,1) - y_history_pid_larger;

% Plot
figure;
subplot(3, 1, 1);
Vtime = (((1:N)*Ts)/60)';
plot(Vtime, y_history_pid_smaller, 'b', 'LineWidth', 1.5);hold on;
plot(Vtime, y_history_pid_normal, 'g', 'LineWidth', 1.5);hold on;
plot(Vtime, y_history_pid_larger, 'c', 'LineWidth', 1.5);hold on;
plot(Vtime, r(:, 1), 'r--', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Process Variable');
title('Process Variable for PI Controller');
legend('Smaller', 'Normal', 'Larger', 'Ref');

subplot(3, 1, 2);
plot(Vtime, y_history_mpc_smaller, 'b', 'LineWidth', 1.5);hold on;
plot(Vtime, y_history_mpc_normal, 'g', 'LineWidth', 1.5);hold on;
plot(Vtime, y_history_mpc_larger, 'c', 'LineWidth', 1.5);hold on;
plot(Vtime, r(:, 1), 'r--', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Process Variable');
title('Process Variable for MPC Controller');
legend('Smaller', 'Normal', 'Larger', 'Ref');

subplot(3, 1, 3);
plot(Vtime, y_history_pid_smaller - y_history_mpc_smaller, 'b', 'LineWidth', 1.5);hold on;
plot(Vtime, y_history_pid_normal - y_history_mpc_normal, 'g', 'LineWidth', 1.5);hold on;
plot(Vtime, y_history_pid_larger - y_history_mpc_larger, 'c', 'LineWidth', 1.5);hold on;
%plot(Vtime, u_history_pid_smaller, 'b', 'LineWidth', 1.5);hold on;
%plot(Vtime, u_history_pid_normal, 'g', 'LineWidth', 1.5);hold on;
%plot(Vtime, u_history_pid_larger, 'c', 'LineWidth', 1.5);hold on;
xlabel('Time'); ylabel('De');
title('Simulation Difference (Accumulated)');
legend('Smaller', 'Normal', 'Larger', 'Ref');

sgtitle('PI and MPC. Different Simulations. Same Ref Signal and Disturbance.')
 

%% Save data as txt file
T_cell = table(Vtime,r(:, 1),y_history_pid_smaller,y_history_mpc_smaller,y_history_pid_smaller - y_history_mpc_smaller, ...
    y_history_pid_normal,y_history_mpc_normal,y_history_pid_normal - y_history_mpc_normal, ...
    y_history_pid_larger,y_history_mpc_larger,y_history_pid_larger - y_history_mpc_larger, ...
    'VariableNames',["time","ref","pi_smaller","mpc_smaller","diff_smaller", ...
                    "pi_normal","mpc_normal","diff_normal", ...
                    "pi_larger","mpc_larger","diff_larger"]);
writetable(T_cell,'txtData/EsimModelError.txt');






