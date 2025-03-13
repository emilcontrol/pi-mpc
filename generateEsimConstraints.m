% This script is made to handle the "normal" simulation but with constraints
% on the process variable.

% Note that for the MPC to operate like in the paper, you will have to
% update the disturbance model in readParamsToWorkspace.m such that 
% d0_ctrl = d0, which means that the disturbance is measured and estimated
% in the MPC to be kept constant.

clear all
[ad,bd,A,B_u,B_d,R,Ts,N,Ndist,K,Ti,Kp,Ki,K_fb,distVec,distVec_ctrl,x0] = readParamsToWorkspace();

% PI controller
PIcontroller = cPI();

[Q, Qf] = calculateQdiscForPI(ad,bd,Ki,Kp,R);

controlHorizon = 60; %200
predHorizon = controlHorizon;

% Reference only for plotting. Always zero.
r = [zeros(N, 1), zeros(N, 1)];

% Initial state
% Minus since we here have y and not e
y_pid = 0;
x_mpc = [0;0];

% History
u_history_pid = zeros(N, 1);
y_history_pid = zeros(N, 1);

u_history_mpc = zeros(N, 1);
x_history_mpc = zeros(N, 2);

% Optimization matrices H and fnx. G and F we need for constraint matrices
[H,fnx_u,fnx_d, G_u, G_d,F] = CalculateQPMtx(A,B_u,B_d,Q,R,Qf,controlHorizon,predHorizon);

% Hard constraints only for MPC. Note that this is constraints on e=-y. Not y.
Xbounds = [-inf  inf ; % x1min, x1max
           -0.5  inf]; % x2min, x2max
Ubounds = [-inf  inf]; % umin, umax

% Constraint matrices.
[Cbar, cnx, Fext, Gdext] = CalculateConstraintMtx(A, B_u, ...
                            Xbounds, Ubounds, controlHorizon, predHorizon, G_u, G_d, F);

% Suppress outputs from quadprog
OptOptions = optimoptions('quadprog', 'Display', 'off','OptimalityTolerance',eps, ...
    'MaxIterations',1000);
warning('off', 'all');

% Simulating
d_ctrl=0;
for k = 1:N
    %Measurable disturbances
    %D = distVec_ctrl(k:k+predHorizon-1);
    D = ones(predHorizon,1)*d_ctrl;

    % control signal for PID
    [u_pid,PIcontroller] = PIcontroller.PI(r(k,2),y_pid,Kp,Ki);

    % Solve QP to find optimal U
    % min U'*H*U/2 + x'*fnx*U subject to Cbar*U <= cnx - Fext*x
    tic
    if k == 1
        U_mpc_old = zeros(controlHorizon,1);
    else
        U_mpc_old = circshift(U_mpc,-1);
    end
    U_mpc = SolveMPC(x_mpc,H,fnx_u,fnx_d, Cbar,cnx,Fext,Gdext,D,U_mpc_old, OptOptions);
    OptTime = toc;

    disp(['Iteration: ', num2str(k), '. Optimization time: ', ...
        num2str(round(OptTime*1000, 1 - floor(log10(abs(OptTime*1000))))), ' ms']);

    u_history_mpc(k) = U_mpc(1,:);
    x_history_mpc(k, :) = x_mpc;

    y_history_pid(k) = y_pid;
    u_history_pid(k) = u_pid;

    %if k > Ndist
    %    d=d0;
    %else
    %    d=0;
    %end
    d=distVec(k);

    % known and measured disturbance vector. We assume all future dist will
    % be like this in next iteration.
    d_ctrl = distVec_ctrl(k);
    
    % Simulate the plant for one step for PID
    y_pid = ad * y_pid + bd * u_pid - B_d(2)*d; 

    % Simulate the plant for one step for MPC
    x_mpc = A * x_mpc + B_u * U_mpc(1,:) + B_d*d; 
end

% Translate MPC x=e back to y=r-e
%y_history_mpc = r(:,1) - x_history_mpc(:, 2);

y_history_mpc = -1*x_history_mpc(:, 2);
%e_history_pid = r(:,1) - y_history_pid(:);

figure;
subplot(2, 1, 1);
Vtime = (((1:N)*Ts)/60)';
plot(Vtime, y_history_pid, 'b', 'LineWidth', 1.5);hold on;
plot(Vtime, y_history_mpc, 'g', 'LineWidth', 1.5);hold on;
plot(Vtime, r(:, 1), 'r--', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Process Variable');
title('Process Variable');
legend('PI', 'MPC');

subplot(2, 1, 2);
plot(Vtime, u_history_pid, 'b', 'LineWidth', 1.5);hold on;
plot(Vtime, u_history_mpc, 'g', 'LineWidth', 1.5);
xlabel('Time'); ylabel('U');
title('Control Signals');
legend('PI', 'MPC');

sgtitle('PI and MPC. Different Simulations. Same Ref Signal and Disturbance.')

%% 
% Save data as txt file
T_cell = table(Vtime,r(:, 1),y_history_pid,y_history_mpc,u_history_pid, ...
    u_history_mpc,'VariableNames',["time","ref","y_pi","y_mpc","u_pi","u_mpc"]);
writetable(T_cell,'txtData/EsimConstraints.txt');
